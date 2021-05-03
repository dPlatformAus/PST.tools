struct Geom_Eigen_Data
{
    std::vector<double> vect; // the eigenvector 0=x, 1=y, 2=z
    double val;	// the eigenvalue
    Geom_Eigen_Data():vect(3,0){};
};
//***********************************************************************************************************************************************************************************
class Geom_Atom : public Atom
{
    public:
		std::vector<double> com_coord; // the centre of mass coordinates, 0=x, 1=y, 2=z
		std::vector<double> pas_coord; // the principal axes system coordinates, 0=x, 1=y, 2=z
		Geom_Atom():com_coord(3,0),pas_coord(3,0){};
};
//***********************************************************************************************************************************************************************************
class Geom_Molecule
{
    public:
        std::string name;
        std::vector<Geom_Atom> atoms; //holds Geom_Atom objects for each atom in the molecule
        double molecular_mass; //sum of atomic masses
        std::vector<double> centre_of_mass_coord; // the centre of mass coordinates, 0=x, 1=y, 2=z
        TNT::Matrix <double> distance_matrix;
        TNT::Array2D <double> inertia_tensor;
        std::vector<Geom_Eigen_Data> inertia_eigendata; //will point to an array of 3 Eigen_Data structs containing the principal moments of inertia (val) and their corresponding axes (vect)
        std::vector<double> inertia_xyz;
        double geometric_mean_I_tensor, arithmetic_mean_I_tensor, bcmplx_tensor, geometric_mean_I_xyz, arithmetic_mean_I_xyz, bcmplx_xyz; //(IA*IB*IC)^(1/3) is the geometric mean
        double geometric_B_tensor, geometric_B_xyz;
        double hindered_rotational_constant, avg_hinderance_rad, avg_hinderance_deg;
        bool no_torsion_tensor, no_2d_rotor_tensor;
        bool no_torsion_xyz, no_2d_rotor_xyz;

        Geom_Molecule(const std::string in_name);
        ~Geom_Molecule(){};
        void add_atom(const std::string &atom_symbol, const double &in_x, const double &in_y, const double &in_z);
        void write_cartesian(FILE *the_file);
        void get_mass();
        double get_sub_mass(const std::vector<unsigned> &target_atoms);
        void get_centre_of_mass_coordinates();
        std::vector<double> get_sub_centre_of_mass_coordinates(const std::vector<unsigned> &target_atoms);
        void get_distance_matrix(FILE *output_file, const bool &output_to_file);
        void write_distance_matrix(FILE *the_file);
        void get_inertia_xyz(); //get moments of inertia about x, y and z coordinates
        void write_inertia_xyz(FILE *the_file);
        void get_inertia_tensor(); //get principal moments of inertia (may not be about x, y and z)
        std::vector<Geom_Eigen_Data> get_sub_inertia_tensor(const std::vector<unsigned> &target_atoms, const std::vector<double> &centre_of_mass); //get principal moments of inertia (may not be about x, y and z) for a subset of atoms
        double get_sub_inertia_axis_vector(const std::vector<unsigned> &target_atoms, const std::vector<double> &centre_of_mass, const std::vector<double> &axis_vector);
        void write_inertia_tensor(FILE *the_file);
        void get_principal_axes_system_coordinates();
        void write_principal_axes_system_coordinates(FILE *the_file);
        void write_xyz_moments_of_inertia(FILE *the_file);
        void write_principal_moments_of_inertia(FILE *the_file);
        void get_simple_properties(double &out_mass, double &out_A, double &out_B, double &out_C);
        double get_angle(const unsigned &a, const unsigned &b, const unsigned &c);
        double get_angle(const unsigned &a, const unsigned &b, const std::vector<double> &point);
        std::vector<double> get_vector(const unsigned &a, const unsigned &b);
        std::vector<double> get_vector(const unsigned &a, const std::vector<double> &point);
        double get_distance_from_com(const unsigned &target_atom);
        double get_distance_from_point(const unsigned &target_atom, const std::vector<double> &point);
        double get_TS_inertia(); //this will have to get smarter, as we may have to get it around an axis perpendicular to the breaking bond or somthing???
};
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Geom_Molecule::Geom_Molecule(const std::string in_name):centre_of_mass_coord(3,0),inertia_tensor(3, 3),inertia_eigendata(3),inertia_xyz(3,0){
    name=in_name;
    no_torsion_tensor=no_2d_rotor_tensor=no_torsion_xyz=no_2d_rotor_xyz=0;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::add_atom(const std::string &atom_symbol, const double &in_x, const double &in_y, const double &in_z){
    Geom_Atom new_atom;
    new_atom.set_atom(atom_symbol,in_x,in_y,in_z);
    atoms.push_back(new_atom);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Geom_Molecule::get_angle(const unsigned &a, const unsigned &b, const unsigned &c) { //using the cosine rule, returns the angle abc defined by 3 atoms where ab and bc define the rays and b is at the vertex, in radians
    double ab, bc, ac, theta;
    ab=sqrt(pow((atoms[a].coord[0]-atoms[b].coord[0]),2)+pow((atoms[a].coord[1]-atoms[b].coord[1]),2)+pow((atoms[a].coord[2]-atoms[b].coord[2]),2));
    bc=sqrt(pow((atoms[b].coord[0]-atoms[c].coord[0]),2)+pow((atoms[b].coord[1]-atoms[c].coord[1]),2)+pow((atoms[b].coord[2]-atoms[c].coord[2]),2));
    ac=sqrt(pow((atoms[a].coord[0]-atoms[c].coord[0]),2)+pow((atoms[a].coord[1]-atoms[c].coord[1]),2)+pow((atoms[a].coord[2]-atoms[c].coord[2]),2));
    theta=acos((ab*ab+bc*bc-ac*ac)/(2*ab*bc));
    return theta;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Geom_Molecule::get_angle(const unsigned &a, const unsigned &b, const std::vector<double> &point) { //using the cosine rule, returns the angle a-b-point defined by 2 atoms and a point in space where a-b and b-point define the rays and b is at the vertex, in radians
    double ab, bc, ac, theta;
    ab=sqrt(pow((atoms[a].coord[0]-atoms[b].coord[0]),2)+pow((atoms[a].coord[1]-atoms[b].coord[1]),2)+pow((atoms[a].coord[2]-atoms[b].coord[2]),2));
    bc=sqrt(pow((atoms[b].coord[0]-point[0]),2)+pow((atoms[b].coord[1]-point[1]),2)+pow((atoms[b].coord[2]-point[2]),2));
    ac=sqrt(pow((atoms[a].coord[0]-point[0]),2)+pow((atoms[a].coord[1]-point[1]),2)+pow((atoms[a].coord[2]-point[2]),2));
    theta=acos((ab*ab+bc*bc-ac*ac)/(2*ab*bc));
    return theta;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<double> Geom_Molecule::get_vector(const unsigned &a, const unsigned &b){ //returns the vector from atom a to atom b
    unsigned i;
    std::vector<double> result;
    for (i=0;i<3;i++) result.push_back(atoms[b].coord[i]-atoms[a].coord[i]);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<double> Geom_Molecule::get_vector(const unsigned &a, const std::vector<double> &point){ //returns the vector from atom a to point, a point in space
    unsigned i;
    std::vector<double> result;
    for (i=0;i<3;i++) result.push_back(point[i]-atoms[a].coord[i]);
    return result;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Geom_Molecule::get_distance_from_com(const unsigned &target_atom){
    unsigned i;
    std::vector<double> atom_to_com;
    get_centre_of_mass_coordinates();//calls get_mass
    for (i=0;i<3;i++) atom_to_com.push_back(atoms[target_atom].coord[i]-centre_of_mass_coord[i]);
    return vector_magnitude(atom_to_com);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Geom_Molecule::get_distance_from_point(const unsigned &target_atom, const std::vector<double> &point){
    unsigned i;
    std::vector<double> atom_to_point;
    for (i=0;i<3;i++) atom_to_point.push_back(atoms[target_atom].coord[i]-point[i]);
    return vector_magnitude(atom_to_point);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Geom_Molecule:: get_TS_inertia(){
    get_centre_of_mass_coordinates();//calls get_mass
    get_inertia_tensor();
    return inertia_eigendata[1].val; //this is only ok for CO or other linear molecules ... more work required for other shapes!
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::get_simple_properties(double &out_mass, double &out_A, double &out_B, double &out_C) {
    get_centre_of_mass_coordinates();//calls get_mass
    get_inertia_tensor();
    out_mass=molecular_mass;
    if (inertia_eigendata[2].val > 0) out_A=I_AMU_SQANGSTROM_TO_B_WAVENUMBER/inertia_eigendata[2].val;
    else out_A=0; //linear molecules will have a 0 moment of inertia along their axis ... this would result in an infinite rotational constant ... linear molecules will just use B
    out_B=I_AMU_SQANGSTROM_TO_B_WAVENUMBER/inertia_eigendata[1].val;
    out_C=I_AMU_SQANGSTROM_TO_B_WAVENUMBER/inertia_eigendata[0].val;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::write_cartesian(FILE *the_file) {
    unsigned i,j;
	fprintf(the_file, "\n%s CARTESIAN COORDINATES:\n", name.c_str());
	fprintf(the_file, "ATOM     MASS     RADIUS(default)   X             Y            Z\n");
    for (i=0;i<atoms.size();i++){
        fprintf(the_file, "%d       %lf   %lf", i, atoms[i].mass, atoms[i].radius);
        for (j=0;j<3;j++) fprintf(the_file, "       %lf", atoms[i].coord[j]);
        fprintf(the_file, "\n");
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::get_distance_matrix(FILE *output_file, const bool &output_to_file) {
    unsigned i,j;
    distance_matrix.newsize(atoms.size()-1, atoms.size()-1);

    for (i=0;i<atoms.size()-1;i++){
        for (j=i+1;j<atoms.size()-1;j++) distance_matrix[j][i]=0;
    }
    for (i=1;i<atoms.size();i++){
        for (j=0;j<i;j++){
            distance_matrix[j][i-1]=sqrt(pow((atoms[i].coord[0]-atoms[j].coord[0]),2)+pow((atoms[i].coord[1]-atoms[j].coord[1]),2)+pow((atoms[i].coord[2]-atoms[j].coord[2]),2));
            if (distance_matrix[j][i-1] <= 0.95){
                if (output_to_file) fprintf(output_file, "Warning: atoms %d and %d are very close!\n", i+1, j+1);
                //because this is a warning, se will always output it to the screen
                fprintf(stderr, "Warning: atoms %d and %d are very close!\n", i+1, j+1);
            }
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::write_distance_matrix(FILE *the_file) {
    unsigned i,j;
    fprintf(the_file, "\nDISTANCES:\n");
    for (i=1;i<atoms.size();i++) fprintf(the_file, "     %d      ", i);
    fprintf(the_file, "\n");
    for (i=0;i<(unsigned)distance_matrix.num_rows();i++){
        fprintf(the_file, "%d    ", i+2);
        for (j=0;j<(unsigned)distance_matrix.num_cols();j++){
            if (distance_matrix[j][i]) fprintf(the_file, "%lf    ", distance_matrix[j][i]);
        }
        fprintf(the_file, "\n");
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::get_principal_axes_system_coordinates() {
    unsigned i,j,k;
    for (i=0;i<atoms.size();i++){
        for (j=0;j<3;j++) {             //j loops through x,y,z of principal axes system coordinates
            atoms[i].pas_coord[j]=0;    //k loops through x,y,z of centre of mass coordinates
            for (k=0;k<3;k++) atoms[i].pas_coord[j]+=atoms[i].com_coord[k]*inertia_eigendata[j].vect[k];
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::write_xyz_moments_of_inertia(FILE *the_file) {
	fprintf(the_file, "\n%s XYZ MOMENTS OF INERTIA (AMU ANGSTROM^2):\n", name.c_str());
	fprintf(the_file, "%lf       %lf       %lf\n", inertia_xyz[0], inertia_xyz[1], inertia_xyz[2]);
	fprintf(the_file, "\nI geometric mean = (IA*IB*IC)^(1/3) = %lf\n", geometric_mean_I_xyz);
	fprintf(the_file, "\nI arithmetic mean = (IA+IB+IC)/3 = %lf\n", arithmetic_mean_I_xyz);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::write_principal_moments_of_inertia(FILE *the_file) {
	fprintf(the_file, "\n%s PRINCIPAL MOMENTS OF INERTIA (AMU ANGSTROM^2):\n", name.c_str());
	fprintf(the_file, "%lf       %lf       %lf\n", inertia_eigendata[0].val, inertia_eigendata[1].val, inertia_eigendata[2].val);
	fprintf(the_file, "\nI geometric mean = (IA*IB*IC)^(1/3) = %lf\n", geometric_mean_I_tensor);
	fprintf(the_file, "\nI arithmetic mean = (IA+IB+IC)/3 = %lf\n", arithmetic_mean_I_tensor);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::write_principal_axes_system_coordinates(FILE *the_file) {
    unsigned i,j;
	fprintf(the_file, "\nPRINCIPAL AXES SYSTEM COORDINATES:\n");
	fprintf(the_file, "ATOM     MASS             X             Y            Z\n");
    for (i=0;i<atoms.size();i++){
        fprintf(the_file, "%d       %lf", i, atoms[i].mass);
        for (j=0;j<3;j++) fprintf(the_file, "       %lf", atoms[i].pas_coord[j]);
        fprintf(the_file, "\n");
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::get_mass() {
    unsigned i;
    molecular_mass=0;
    for (i=0;i<atoms.size();i++) molecular_mass+=atoms[i].mass;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Geom_Molecule::get_sub_mass(const std::vector<unsigned> &target_atoms) {
    unsigned i;
    double sub_mass=0;
    for (i=0;i<target_atoms.size();i++) sub_mass+=atoms[ target_atoms[i] ].mass;
    return sub_mass;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::get_centre_of_mass_coordinates() {
    unsigned i,j;
    std::vector<double> mass_component(3,0);
    get_mass();
    for (i=0;i<atoms.size();i++) for (j=0;j<3;j++) mass_component[j]+=(atoms[i].mass*atoms[i].coord[j]);
    for (j=0;j<3;j++) centre_of_mass_coord[j]=mass_component[j]/molecular_mass;
    for (i=0;i<atoms.size();i++){
        for (j=0;j<3;j++) atoms[i].com_coord[j]=atoms[i].coord[j]-centre_of_mass_coord[j];
    }
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<double> Geom_Molecule::get_sub_centre_of_mass_coordinates(const std::vector<unsigned> &target_atoms) {
    unsigned i,j;
    std::vector<double> mass_component(3,0);
    std::vector<double> sub_centre_of_mass_coord(3,0);
    double sub_mass=get_sub_mass(target_atoms);
    for (i=0;i<target_atoms.size();i++) for (j=0;j<3;j++) mass_component[j] += ( atoms[ target_atoms[i] ].mass * atoms[ target_atoms[i] ].coord[j] );
    for (j=0;j<3;j++) sub_centre_of_mass_coord[j]=mass_component[j]/sub_mass;
    return sub_centre_of_mass_coord;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::get_inertia_xyz() {
    unsigned i;
    //double ordercheck1, ordercheck2, temp_inertia;
    for (i=0;i<3;i++) inertia_xyz[i]=0; //reset moments of inertia
    for (i=0;i<atoms.size();i++){
        inertia_xyz[0]+=atoms[i].mass*(pow(atoms[i].com_coord[1],2)+pow(atoms[i].com_coord[2],2)); //Ixx+=mass*(y^2+z^2)
        inertia_xyz[1]+=atoms[i].mass*(pow(atoms[i].com_coord[0],2)+pow(atoms[i].com_coord[2],2)); //Iyy+=mass*(x^2+z^2)
        inertia_xyz[2]+=atoms[i].mass*(pow(atoms[i].com_coord[0],2)+pow(atoms[i].com_coord[1],2)); //Izz+=mass*(x^2+y^2)
    }
    /* ************* we don't want to do this because the Z inertia is the torsion!! *********
    std::sort(inertia_xyz, inertia_xyz+3);
    //We sorted ascending for the (ordercheck2<ordercheck1) swap to work
    ordercheck1=inertia_xyz[1]/inertia_xyz[0];
    ordercheck2=inertia_xyz[2]/inertia_xyz[1];
    // ************************* SORT SO THAT THE LEAST SIMILAR ONE IS ISOLATED in [2]! *************************
    if (ordercheck2<ordercheck1){
        temp_inertia=inertia_xyz[2];
        inertia_xyz[2]=inertia_xyz[0];
        inertia_xyz[0]=temp_inertia;
    }*/
    geometric_mean_I_xyz=cbrt(inertia_xyz[0]*inertia_xyz[1]*inertia_xyz[2]); //(Ix*Iy*Iz)^(1/3)
    arithmetic_mean_I_xyz=(inertia_xyz[0]+inertia_xyz[1]+inertia_xyz[2])/3; //(Ix+Iy+Iz)/3
    bcmplx_xyz=cbrt((16.86/inertia_xyz[0]) * (16.86/inertia_xyz[1]) * (16.86/inertia_xyz[2]) ); //(16.86/Ix * 16.86/Iy * 16.86/Iz)^(1/3)
    if (inertia_xyz[2] <= 0.05) no_torsion_xyz=1;
    if (inertia_xyz[0] <= 0.05 || inertia_xyz[1] <= 0.05) no_2d_rotor_xyz=1;
    //I_AMU_SQANGSTROM_TO_B_WAVENUMBER is a conversion factor between AMU.A^2 to CM^-1 ... turning I into B  ...
    geometric_B_xyz=I_AMU_SQANGSTROM_TO_B_WAVENUMBER/sqrt(inertia_xyz[0]*inertia_xyz[1]);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::get_inertia_tensor() {
    /*---------------------------
    The inertia tensor looks like this
        | Ixx  Ixy  Ixz |
        | Iyx  Iyy  Iyz |
        | Izx  Izy  Izz |
    ---------------------------*/
    unsigned i,j;
    TNT::Array2D <double> inertia_eigenvectors(3,3);
    TNT::Array1D <double> inertia_eigenvalues(3);
    double Ixx, Iyy, Izz, Ixy, Ixz, Iyz, eigencheck1, eigencheck2;
    Geom_Eigen_Data temp_eigen_data;
    Ixx=Iyy=Izz=Ixy=Ixz=Iyz=0;
    for (i=0;i<atoms.size();i++){
        Ixx+=atoms[i].mass*(pow(atoms[i].com_coord[1],2)+pow(atoms[i].com_coord[2],2)); //Ixx+=mass*(y^2+z^2)
        Iyy+=atoms[i].mass*(pow(atoms[i].com_coord[0],2)+pow(atoms[i].com_coord[2],2)); //Iyy+=mass*(x^2+z^2)
        Izz+=atoms[i].mass*(pow(atoms[i].com_coord[0],2)+pow(atoms[i].com_coord[1],2)); //Izz+=mass*(x^2+y^2)
        Ixy-=atoms[i].mass*atoms[i].com_coord[0]*atoms[i].com_coord[1]; //Ixy-=mass*x*y
        Ixz-=atoms[i].mass*atoms[i].com_coord[0]*atoms[i].com_coord[2]; //Ixz-=mass*x*z
        Iyz-=atoms[i].mass*atoms[i].com_coord[1]*atoms[i].com_coord[2]; //Iyz-=mass*y*z
    }
    inertia_tensor[0][0]=Ixx;
    inertia_tensor[1][1]=Iyy;
    inertia_tensor[2][2]=Izz;
    inertia_tensor[0][1]=inertia_tensor[1][0]=Ixy;
    inertia_tensor[0][2]=inertia_tensor[2][0]=Ixz;
    inertia_tensor[1][2]=inertia_tensor[2][1]=Iyz;

    JAMA::Eigenvalue <double> eigens(inertia_tensor);
    eigens.getV(inertia_eigenvectors); //results are slightly different from old code, MJTJ said that a factor of -1 difference is common as it is algorith dependent, so it will be consistent, and we may need to modify the results to match ...
    eigens.getRealEigenvalues(inertia_eigenvalues);

    for (i=0;i<3;i++){
        inertia_eigendata[i].val=abs(inertia_eigenvalues[i]);
        for (j=0;j<3;j++) inertia_eigendata[i].vect[j]=inertia_eigenvectors[j][i];
    }
    //JAMA::Eigenvalue returns the results in ascending order, which is what we need for the (eigencheck2<eigencheck1) swap to work
    eigencheck1=inertia_eigendata[1].val/inertia_eigendata[0].val;
    eigencheck2=inertia_eigendata[2].val/inertia_eigendata[1].val;
    //************************* SORT SO THAT THE LEAST SIMILAR ONE IS ISOLATED in [2]! *************************
    if (eigencheck2<eigencheck1){
        temp_eigen_data=inertia_eigendata[2];
        inertia_eigendata[2]=inertia_eigendata[0];
        inertia_eigendata[0]=temp_eigen_data;
    }
    arithmetic_mean_I_tensor=(inertia_eigendata[0].val+inertia_eigendata[1].val+inertia_eigendata[2].val)/3; //  (IA+IB+IC)/3
    geometric_mean_I_tensor=cbrt(inertia_eigendata[0].val*inertia_eigendata[1].val*inertia_eigendata[2].val); //   (IA*IB*IC)^(1/3)
    bcmplx_tensor=cbrt( (16.86/inertia_eigendata[0].val) * (16.86/inertia_eigendata[1].val) * (16.86/inertia_eigendata[2].val) ); //(16.86/IA * 16.86/IB * 16.86 /IC)^(1/3)
    if (inertia_eigendata[2].val <= 0.05) no_torsion_tensor=1;
    if (inertia_eigendata[0].val <= 0.05 || inertia_eigendata[1].val <= 0.05) no_2d_rotor_tensor=1;
    //I_AMU_SQANGSTROM_TO_B_WAVENUMBER is a conversion factor between AMU.A^2 to CM^-1 ... turns I into B
    geometric_B_tensor=I_AMU_SQANGSTROM_TO_B_WAVENUMBER/sqrt(inertia_eigendata[0].val*inertia_eigendata[1].val);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<Geom_Eigen_Data> Geom_Molecule::get_sub_inertia_tensor(const std::vector<unsigned> &target_atoms, const std::vector<double> &centre_of_mass){
    /*---------------------------
    The inertia tensor looks like this
        | Ixx  Ixy  Ixz |
        | Iyx  Iyy  Iyz |
        | Izx  Izy  Izz |
    ---------------------------*/
    unsigned i,j;
    TNT::Array2D <double> sub_inertia_tensor(3, 3);
    TNT::Array2D <double> inertia_eigenvectors(3,3);
    TNT::Array1D <double> inertia_eigenvalues(3);
    double Ixx, Iyy, Izz, Ixy, Ixz, Iyz, eigencheck1, eigencheck2;
    Geom_Eigen_Data temp_eigen_data;
    std::vector<Geom_Eigen_Data> inertia_results(3); //will point to an array of 3 Eigen_Data structs containing the principal moments of inertia (val) and their corresponding axes (vect)
    Ixx=Iyy=Izz=Ixy=Ixz=Iyz=0;
    for (i=0;i<target_atoms.size();i++){
        Ixx+=atoms[target_atoms[i]].mass*(pow(atoms[target_atoms[i]].coord[1]-centre_of_mass[1],2)+pow(atoms[target_atoms[i]].coord[2]-centre_of_mass[2],2)); //Ixx+=mass*(y^2+z^2)
        Iyy+=atoms[target_atoms[i]].mass*(pow(atoms[target_atoms[i]].coord[0]-centre_of_mass[0],2)+pow(atoms[target_atoms[i]].coord[2]-centre_of_mass[2],2)); //Iyy+=mass*(x^2+z^2)
        Izz+=atoms[target_atoms[i]].mass*(pow(atoms[target_atoms[i]].coord[0]-centre_of_mass[0],2)+pow(atoms[target_atoms[i]].coord[1]-centre_of_mass[1],2)); //Izz+=mass*(x^2+y^2)
        Ixy-=atoms[target_atoms[i]].mass*(atoms[target_atoms[i]].coord[0]-centre_of_mass[0])*(atoms[target_atoms[i]].coord[1]-centre_of_mass[1]); //Ixy-=mass*x*y
        Ixz-=atoms[target_atoms[i]].mass*(atoms[target_atoms[i]].coord[0]-centre_of_mass[0])*(atoms[target_atoms[i]].coord[2]-centre_of_mass[2]); //Ixz-=mass*x*z
        Iyz-=atoms[target_atoms[i]].mass*(atoms[target_atoms[i]].coord[1]-centre_of_mass[1])*(atoms[target_atoms[i]].coord[2]-centre_of_mass[2]); //Iyz-=mass*y*z
    }
    sub_inertia_tensor[0][0]=Ixx;
    sub_inertia_tensor[1][1]=Iyy;
    sub_inertia_tensor[2][2]=Izz;
    sub_inertia_tensor[0][1]=sub_inertia_tensor[1][0]=Ixy;
    sub_inertia_tensor[0][2]=sub_inertia_tensor[2][0]=Ixz;
    sub_inertia_tensor[1][2]=sub_inertia_tensor[2][1]=Iyz;
    JAMA::Eigenvalue <double> eigens(sub_inertia_tensor);
    eigens.getV(inertia_eigenvectors); //results are slightly different from old code, MJTJ said that a factor of -1 difference is common as it is algorith dependent, so it will be consistent, and we may need to modify the results to match ...
    eigens.getRealEigenvalues(inertia_eigenvalues);
    for (i=0;i<3;i++){
        inertia_results[i].val=abs(inertia_eigenvalues[i]);
        for (j=0;j<3;j++) inertia_results[i].vect[j]=inertia_eigenvectors[j][i];
    }
    //JAMA::Eigenvalue returns the results in ascending order, which is what we need for the (eigencheck2<eigencheck1) swap to work
    eigencheck1=inertia_results[1].val/inertia_results[0].val;
    eigencheck2=inertia_results[2].val/inertia_results[1].val;
    //************************* SORT SO THAT THE LEAST SIMILAR ONE IS ISOLATED in [2]! *************************
    if (eigencheck2<eigencheck1){
        temp_eigen_data=inertia_results[2];
        inertia_results[2]=inertia_results[0];
        inertia_results[0]=temp_eigen_data;
    }
    return inertia_results;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Geom_Molecule::get_sub_inertia_axis_vector(const std::vector<unsigned> &target_atoms, const std::vector<double> &centre_of_mass, const std::vector<double> &axis_vector){
    unsigned i;
    double angle, distance, perpendicular_distance, I_com_dir(0);
    std::vector<double> com_vector;
    for (i=0;i<target_atoms.size();i++) {
        com_vector=get_vector(target_atoms[i], centre_of_mass);
        distance=vector_magnitude(com_vector);
        angle=vector_angle(com_vector, axis_vector);
        perpendicular_distance=distance*sin(angle);
        I_com_dir+=atoms[target_atoms[i]].mass*pow(perpendicular_distance,2);
    }
    return I_com_dir;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::write_inertia_tensor(FILE *the_file) {
    fprintf(the_file, "\nFOR: %s\n", name.c_str());
    fprintf(the_file, "\nINERTIA TENSOR: \n");
    fprintf(the_file, "%lf %lf %lf \n", inertia_tensor[0][0], inertia_tensor[0][1], inertia_tensor[0][2]);
    fprintf(the_file, "%lf %lf %lf \n", inertia_tensor[1][0], inertia_tensor[1][1], inertia_tensor[1][2]);
    fprintf(the_file, "%lf %lf %lf \n", inertia_tensor[2][0], inertia_tensor[2][1], inertia_tensor[2][2]);
    fprintf(the_file, "\nEIGENVALUES: 1=%lf 2=%lf 3=%lf \n", inertia_eigendata[0].val, inertia_eigendata[1].val, inertia_eigendata[2].val);
    fprintf(the_file, "\nEIGENVECTORS:   x   y   z\n");
    fprintf(the_file, "1  %lf %lf %lf \n", inertia_eigendata[0].vect[0], inertia_eigendata[0].vect[1], inertia_eigendata[0].vect[2]);
    fprintf(the_file, "2  %lf %lf %lf \n", inertia_eigendata[1].vect[0], inertia_eigendata[1].vect[1], inertia_eigendata[1].vect[2]);
    fprintf(the_file, "3  %lf %lf %lf \n", inertia_eigendata[2].vect[0], inertia_eigendata[2].vect[1], inertia_eigendata[2].vect[2]);
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void Geom_Molecule::write_inertia_xyz(FILE *the_file) {
    fprintf(the_file, "\nFOR: %s\n", name.c_str());
    fprintf(the_file, "\nMOMENT OF INERTIA: x=%lf y=%lf z=%lf \n", inertia_xyz[0], inertia_xyz[1], inertia_xyz[2]);
}

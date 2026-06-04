#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

int main(int argc, char *argv[])
{
    double x,y,z,ene,index,force;
    double thea,delx,dely,delz,rxy,r,deltazmc,MC,MO,zmc,deltaz;
    double PI=3.1415926;
    int n,na,natom,ncount,indexn;
    int step,nline,nstep;
    char **element;
    double **fracmat,*fracmat_C,*fracmat_O;
    std::string line;
    std::ifstream ifs;
    std::ofstream ofs;

    MC = 12.00000;
    MO = 15.99491;

    if (argc != 6){
	std::cout<<"Usage: ./exec inputfile(xyzformat) nlineperbolck outfile deltazmc equzmc"<<std::endl;
	exit(1);
    }

    ifs.open(argv[1],std::ios::in);
    nline = atoi(argv[2]);
    ofs.open(argv[3],std::ios::out);
    natom = nline - 2;
    deltazmc = atof(argv[4]);
    zmc = atof(argv[5]);

    fracmat = new double *[natom];
    for (int i = 0;i < natom;i++)
    {
        fracmat[i] = new double [3];
    }

    fracmat_C = new double [3];
    fracmat_O = new double [3];

    element = new char *[natom];
    for (int i = 0;i < natom;i++)
    {
        element[i] = new char [3];
    }
    // std::cout<<"allocate well"<<std::endl;

    if (!ifs.is_open())
    {
        std::cout<<"No such file !"<<std::endl;
        exit(0);
    }

    nstep = 0;
    step = 0;
    ofs<<"step       theta"<<std::endl;
    while (getline(ifs,line))
    {
        ncount = nstep % nline;
        std::stringstream linetemp(line);
        if (ncount == 0)
        {
            linetemp >> na;
        }else if (ncount == 1){
            linetemp >> ene;
        }else{
            linetemp>>element[ncount-2]>>x>>y>>z;
            // std::cout<<element[ncount-2]<<x<<y<<z<<std::endl;
            fracmat[ncount-2][0] = x;
            fracmat[ncount-2][1] = y;
            fracmat[ncount-2][2] = z;
            if (ncount == 34)
            {
                fracmat_C[0] = x;
                fracmat_C[1] = y;
                fracmat_C[2] = z;
                //std::cout<<element[ncount-2]<<x<<y<<z<<std::endl;
            }
            if (ncount == 35)
            {
                fracmat_O[0] = x;
                fracmat_O[1] = y;
                fracmat_O[2] = z;
                //std::cout<<element[ncount-2]<<x<<y<<z<<std::endl;
            }
            if (ncount == (nline - 1) )
            {
                delx = fracmat_O[0] - fracmat_C[0];
                dely = fracmat_O[1] - fracmat_C[1];
                delz = fracmat_O[2] - fracmat_C[2];
                r = sqrt(pow(delx,2)+pow(dely,2)+pow(delz,2));
                rxy = sqrt(pow(delx,2)+pow(dely,2));
                thea = acos(delz/r);
                deltaz = (fracmat_C[2]*MC + fracmat_O[2]*MO)/(MC+MO)-zmc;
                // thea = asin(rxy/r);
                // thea = atan(rxy/delz);
                // std::cout<<std::right<<std::setw(10)<<step<<std::fixed<<std::setprecision(5)
                //          <<std::right<<std::setw(10)<<180.0/PI*thea
                //          <<std::right<<std::setw(10)<<(180.0/PI*thea)/sin(thea)<<std::endl;
                // std::cout<<std::right<<std::setw(10)<<step<<std::fixed<<std::setprecision(5)
                // <<std::right<<std::setw(10)<<delx
                // <<std::right<<std::setw(10)<<dely
                // <<std::right<<std::setw(10)<<delz
                // <<std::right<<std::setw(10)<<rxy<<std::endl;
                if (fabs(deltaz) < deltazmc){
                ofs<<std::right<<std::setw(10)<<step<<std::fixed<<std::setprecision(5)
                    <<std::right<<std::setw(10)<<180.0/PI*thea<<std::endl;
                }
                    // <<std::right<<std::setw(10)<<(180.0/PI*thea)/sin(thea)<<std::endl;
                step++;
            }
        }
        nstep++;
    }
    for (int i = 0;i < natom;i++)
    {
        delete fracmat[i];
    }
    delete fracmat;
    delete fracmat_C;
    delete fracmat_O;
    delete element;
    ifs.close();
    ofs.close();
    
    return 0;
}

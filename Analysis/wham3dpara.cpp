#include <iostream>
#include <string>
#include <cmath>
#include <sstream>
//#include <cerrno>
//#include <cassert>
//#include <cctype>
#include <iomanip>
#include <fstream>
#include <omp.h>
//to compile : g++/icpc wham3dpara.cpp -fopenmp -O2 -o execfile

#define MAXNITER 10000
#define RADIANS   6.28318530717959
#define PI 3.141592653589793
#define KBOLTZ 0.0019872064689 
#define LINESIZE 1024
#define TOL 1.e-4
//#define EQU(a,b) (strncasecmp(a,b,strlen(b))==0)
#define COMMAND_LINE "wham metafile freefile temp minx maxx num_binsx miny maxy num_binsy minz maxz num_binsz PxNo|Px=xx PyNo|Py=xx PzNo|Pz=xx nproc\n"

class wham_info
{
    public:
        int periodic[3];
        double period[3];
        double hist_min[3],hist_max[3],bin_width[3];
        int num_bins[3];
        int num_windows;
        double T;
        double (*bias_locations)[3];
        double (*spring_constants)[3];
        double *F;
        double *kT;
        double *partition;
        double ***X,***rewX;
        int *nt;
        double *ebf,*ebfo,*ebf2,*fact,***ebw;
        //char metaFileName[LINESIZE],fepFileName[LINESIZE],(*dataFileName)[LINESIZE];
        char *metaFileName,*fepFileName;
        std::string *dataFileName;
        double ***prob,***pmf;
};

// class pmf_rew_data
// {
//     public:
//         int num_bins[3];
//         double hist_min[3], hist_max[3], bin_width[3];
//         //int nDim;
//         //char metaFileName[LINESIZE], whamFileName[LINESIZE];
//         // FILE *metaFile, *whamFile;
//         char *metaFile,*whamFile;
//         int num_windows;
//         double (*bias_pos)[3], (*force_const)[3];
//         double *part;
//         double *pmf1, *prob1, **histogram1, **histogramb1;
//         double **pmf2, **prob2, ***histogram2, ***histogramb2;
//         double ***pmf3, ***prob3, ****histogram3, ****histogramb3;
// };

double calc_bias_pot(class wham_info *h, int index, double *coor,int BIAS);
int get_numwindows(char *file);
int is_metadata(char *line);

int main(int argc, char *argv[])
{
    int i, j, k, kj, kjl, l, m, n, m1, m2, m3, n1, n2, n3, i1, i2, i3, j1, j2 , j3, k1, k2, k3,nthread;//20240806
    int current_window, vals, nline, conv;
    char *c;//, filename[LINESIZE], filename2[LINESIZE];
    double loc, loc1, loc2, loc3, spring, spring1, spring2, spring3, correl_time, temp;//20240723 loc3 and spring3
    double beta, e_rstr, coor[3], step, e_sits;//20240723 coor[3]
    double ebfk, bottom, summ, delta;
    double ***histogram, denominator, fmin, fmin1, fmin2, fmin3, fk, nk;
    double a, b, nkValues, betaK, e, x, y, z;
    //double gft[2000], gftt, a0, b0, a1, b1;
    int index, index1, index2, index3;
    class wham_info *wham;
    class pmf_rew_data *pmf;
    std::ifstream file,file2;
    std::ifstream metaFile;
    std::string line,filename;

    if (argc < 3)
    {
        std::cout<<COMMAND_LINE<<std::endl;
    }
    //flush(stdout);

    wham = new class wham_info [1];
    // pmf = new class pmf_rew_data [1];

    //wham->T = 300.;

    
    for(i=0; i<3; i++) {//20240812
        wham->hist_min[i] = -PI;
        wham->hist_max[i] = PI;
        wham->num_bins[i] = 60;
        wham->periodic[i] = 1;
        wham->period[i] = 2*PI;

        // pmf->hist_min[i] = -PI;
        // pmf->hist_max[i] = PI;
        // pmf->num_bins[i] = 60;
    }
    for (i=0;i<argc;i++)
    {
        std::cout<<std::setw(10)<<argv[i];
    }
    std::cout<<std::endl;

    wham->metaFileName=argv[1];
    wham->fepFileName=argv[2];
    wham->T = atof(argv[3]);

    wham->hist_min[0] = atof(argv[4]);
    wham->hist_max[0] = atof(argv[5]);
    wham->num_bins[0] = atoi(argv[6]);

    wham->hist_min[1] = atof(argv[7]);
    wham->hist_max[1] = atof(argv[8]);
    wham->num_bins[1] = atoi(argv[9]);

    wham->hist_min[2] = atof(argv[10]);
    wham->hist_max[2] = atof(argv[11]);
    wham->num_bins[2] = atoi(argv[12]);

    if (argv[13] == "PxNo")
    {
        wham->period[0] = 0;
    }else{
        wham->period[0] = atof(argv[13]);
    }

    if (argv[14] == "PyNo")
    {
        wham->period[0] = 0;
    }else{
        wham->period[0] = atof(argv[14]);
    }

    if (argv[15] == "PzNo")
    {
        wham->period[0] = 0;
    }else{
        wham->period[0] = atof(argv[15]);
    }

    nthread = atoi(argv[16]);

    if (argc > 17)
    {
        std::cout<<"Invalid argument number"<<std::endl;
        exit(0);
    }

    /*calculate bin_width*/
    for(i=0; i<3; i++) {
        wham->bin_width[i] = (wham->hist_max[i] - wham->hist_min[i])/(double)wham->num_bins[i];
    }

    // for(i=0; i<3; i++) {
    //     pmf->bin_width[i] = (pmf->hist_max[i] - pmf->hist_min[i])/(double)pmf->num_bins[i];
    // }

    for(i=0; i<3; i++) {
        if(wham->periodic[i] == 0) wham->period[i] = 1.e200;
    }

    std::cout<<"WHAM 3D"<<std::endl;
    std::cout<<"metaFile: "<<wham->metaFileName<<std::endl;
    std::cout<<"Temperature: "<<wham->T<<std::endl;
    for (i=0;i<3;i++)
    {
        std::cout<<"hist_min hist_max hist_bins --"<<i<<std::endl;
        std::cout<<std::setw(10)<<wham->hist_min[i]<<std::setw(10)<<wham->hist_max[i]<<std::setw(10)<<wham->num_bins[i]<<std::endl;
    }
    //flush(stdout);

    wham->num_windows = get_numwindows(wham->metaFileName);
    std::cout<<"numbers of windows used: "<<wham->num_windows<<std::endl;
    //flush(stdout);

    /*read data from metadafile */
    metaFile.open(wham->metaFileName,std::ios::in);
    if(!metaFile.is_open())
    {
        std::cout<<"No such File !"<<std::endl;
        exit(0);
    }

    wham->F = new double [wham->num_windows];
    wham->kT = new double [wham->num_windows];
    wham->partition = new double [wham->num_windows];
    wham->nt = new int [wham->num_windows];
    wham->dataFileName = new std::string [wham->num_windows];
    wham->bias_locations = new double [wham->num_windows][3];
    wham->spring_constants = new double [wham->num_windows][3];
    wham->prob = new double **[wham->num_bins[0]];
    wham->pmf = new double **[wham->num_bins[0]];
    for (i=0;i<wham->num_bins[0];i++)
    {
        wham->prob[i] = new double *[wham->num_bins[1]];
        wham->pmf[i] = new double *[wham->num_bins[1]];
        for (j=0;j<wham->num_bins[1];j++)
        {
            wham->prob[i][j] = new double [wham->num_bins[2]];
            wham->pmf[i][j] = new double [wham->num_bins[2]];
        }
    }
   
    current_window = 0;
    while (getline(metaFile,line))
    {
        std::stringstream iss(line);
        iss>>filename>>loc1>>loc2>>loc3>>spring1>>spring2>>spring3;//>>correl_time>>temp;
        wham->dataFileName[current_window]=filename;
        wham->bias_locations[current_window][0] = loc1;
        wham->bias_locations[current_window][1] = loc2;
        wham->bias_locations[current_window][2] = loc3;
        wham->spring_constants[current_window][0] = spring1;
        wham->spring_constants[current_window][1] = spring2;
        wham->spring_constants[current_window][2] = spring3;

        current_window++;

    }

    if (current_window != wham->num_windows){
        std::cout<<"Inconsistent no. of windows"<<current_window<<" "<<wham->num_windows<<std::endl;
        exit(0);
    }else{
        std::cout<<"read file finished!"<<std::endl;
    }

    for (i = 0;i<wham->num_windows;i++)
    {
        file.open(wham->dataFileName[i],std::ios::in);
        if (!file.is_open())
        {
            std::cout<<"can not open the file: "<<wham->dataFileName[i]<<std::endl;
            exit(0);
        }
        nline = 0;
        while(getline(file,line)){
            nline++;
        }
        wham->nt[i] = nline;
        wham->partition[i] = nline;
        file.close();
    }
    std::cout<<"count line finished! "<<std::endl;

    wham->ebw = new double **[wham->num_windows];
    for (i=0;i<wham->num_windows;i++)
    {
        wham->ebw[i] = new double *[wham->nt[i]];
        for (j = 0;j<wham->nt[i];j++)
        {
            wham->ebw[i][j] = new double [wham->num_windows];
        }
    }

    wham->X = new double **[wham->num_windows];
    for (i = 0 ;i<wham->num_windows;i++)
    {
        wham->X[i] = new double *[wham->nt[i]];
        for (j=0;j<wham->nt[i];j++)
        {
            wham->X[i][j] = new double [3];
        }
    }

    wham->rewX = new double **[wham->num_windows];
    for (i = 0 ;i<wham->num_windows;i++)
    {
        wham->rewX[i] = new double *[wham->nt[i]];
        for (j=0;j<wham->nt[i];j++)
        {
            wham->rewX[i][j] = new double [3];
        }
    }

    beta = 1./(wham->T*KBOLTZ);
    if (wham->kT[0] > 1.e-3)
    {
        for(i=0;i<wham->num_windows;i++)
        {
            wham->kT[i] = 1./(wham->kT[i]*KBOLTZ);
            std::cout<<"REMD temp: "<<i<<" "<<wham->kT[i]<<std::endl;
        }
    }
    //flush(stdout);

    for (i=0;i<wham->num_windows;i++)
    {
        file.open(wham->dataFileName[i],std::ios::in);
        if(!file.is_open())
        {
            std::cout<<"Can not open the file: "<<wham->dataFileName[i]<<std::endl;
            exit(0);
        }
        l = 0;
        while(getline(file,line))
        {
            std::stringstream isse(line);
            isse>>step>>coor[0]>>coor[1]>>coor[2];
            wham->X[i][l][0] = coor[0];
            wham->X[i][l][1] = coor[1];
            wham->X[i][l][2] = coor[2];
            for (k = 0;k<wham->num_windows;k++)
            {
                e_rstr = calc_bias_pot(wham,k,coor,0);
                wham->ebw[i][l][k] = exp(-beta*e_rstr);
            }
            l++;
        }
        if (l != wham->nt[i])
        {
            std::cout<<"Inconsistent no. of lines in data file: "<<wham->dataFileName[i]<<" "<<wham->nt[i]<<std::endl;
        }
        file.close();    
    }
    std::cout<<"read data finished!"<<std::endl;

    wham->ebf = new double [wham->num_windows];
    wham->ebfo = new double [wham->num_windows];
    wham->fact = new double [wham->num_windows];
    wham->ebf2 = new double [wham->num_windows];
    for(i=0; i<wham->num_windows; i++) {
        wham->F[i] = 0.;
        wham->ebf[i] = exp(beta*wham->F[i]);
        wham->ebfo[i] = wham->ebf[i];
        wham->fact[i] = wham->nt[i]*wham->ebf[i];
    }///
    //wham iteration
    for (n = 0;n<=MAXNITER;n++)
    {
        for (k=0;k<wham->num_windows;k++)
        {
            ebfk = 0.;
            # pragma omp parallel for num_threads(nthread) \
                reduction(+: ebfk) schedule(guided)
            for (i = 0;i<wham->num_windows;i++)
            {
                for (l = 0;l<wham->nt[i];l++)
                {
                    bottom = 0.;
                    for (j = 0; j<wham->num_windows;j++)
                    {
                        bottom += wham->ebw[i][l][j]*wham->fact[j];
                    }
                    ebfk += wham->ebw[i][l][k]/bottom;
                }
            }
            wham->ebf2[k] = ebfk;
            wham->ebf[k] = 1./(wham->ebf[0]*ebfk);
            wham->fact[k] = wham->nt[k]*wham->ebf[k];
        }

        conv=1;
        summ = 0.;
        for(k = 0;k<wham->num_windows;k++)
        {
            wham->ebf[k] = wham->ebf2[0]/wham->ebf2[k];
            delta = fabs(1./beta*log(wham->ebf[k]/wham->ebfo[k]));

            if(delta >= TOL) conv = 0;
            wham->ebfo[k] = wham->ebf[k];
            summ += delta;
        }

        std::cout<<n<<" "<<summ<<std::endl;
        if (n%100 == 0)
        {
            for(i = 0;i<wham->num_windows;i++)
            {
                std::cout<<"$$$$ "<<1./beta*log(wham->ebf[i])<<std::endl;
            }
        }
        if(conv > 0) break;
        //flush(stdout);
    }
    delete wham->fact;
    delete wham->ebfo;
    delete wham->ebf2;


    std::cout<<"converged F[k]:"<<std::endl;
    for (i = 0;i<wham->num_windows;i++)
    {
        wham->F[i] = 1./beta*log(wham->ebf[i]);
        std::cout<<std::right<<std::setw(12)<<wham->nt[i]<<
        std::fixed<<std::setprecision(5)<<std::right<<std::setw(12)<<wham->F[i]<<std::endl;
    }
    //flush(stdout);


    histogram = new double **[wham->num_bins[0]];
    for (i = 0;i<wham->num_bins[0];i++)
    {
        histogram[i] = new double *[wham->num_bins[1]];
        for(j = 0;j<wham->num_bins[1];j++)
        {
            histogram[i][j] = new double [wham->num_bins[2]];
        }
    }

    for (i=0; i<wham->num_windows; i++) {//20240806
        for(k=0; k<wham->num_bins[0]; k++){
           for(kj=0; kj<wham->num_bins[1]; kj++){
                for (int kjl=0;kjl<wham->num_bins[2];kjl++){
                    histogram[k][kj][kjl] = 0.;
                }
            }
        }
        
        n1 = 0;
        for(j = 0;j<wham->nt[i];j++)
        {
            k1 = 0;k2=0;k3=0;
            if(wham->X[i][j][0] >= wham->hist_min[0] && wham->X[i][j][0] < wham->hist_max[0] 
             && wham->X[i][j][1] >= wham->hist_min[1] && wham->X[i][j][1] < wham->hist_max[1] 
             && wham->X[i][j][2] >= wham->hist_min[2] && wham->X[i][j][2] < wham->hist_max[2] ) {
                        index1 = 
                        (int) ((wham->X[i][j][0] - wham->hist_min[0])/wham->bin_width[0]);
                        index2 = 
                        (int) ((wham->X[i][j][1] - wham->hist_min[1])/wham->bin_width[1]);
                        index3 = 
                        (int) ((wham->X[i][j][2] - wham->hist_min[2])/wham->bin_width[2]);
                        k3 = 1;
                        k2 = 1;
                        k1 = 1;
                        n1++;
            }

            denominator = 0.;
                        coor[0] = wham->X[i][j][0];
            coor[1] = wham->X[i][j][1];
            coor[2] = wham->X[i][j][2];//20240806
            for (k=0; k<wham->num_windows; k++) {
                e_rstr = calc_bias_pot(wham, k, coor, 0);
                fk = wham->F[k];
                nk = wham->nt[k];
                denominator += nk*exp(-beta*(e_rstr - fk));
            }
            
            if(k1 == 1 && k2 == 1 && k3 == 1) {//20240806
                histogram[index1][index2][index3] += 1./denominator;
            }
        }
        for(k=0; k<wham->num_bins[0]; k++) {
            for(kj=0; kj<wham->num_bins[1]; kj++) {//20240806
                for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                    wham->prob[k][kj][kjl] += histogram[k][kj][kjl];//20240806
                }
            }
        }
        /*debug*/
        //fclose(file2);
        /*end debug*/
    }

    /*normalization*/
    summ = 0.;
    for(k=0; k<wham->num_bins[0]; k++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                summ += wham->prob[k][kj][kjl];
            }//20240806
        }
    }

    for(k=0; k<wham->num_bins[0]; k++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                wham->prob[k][kj][kjl] /= summ;
            }//20240806
        }
    }

    for(k=0; k<wham->num_bins[0]; k++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                std::cout<<"prob: "<<std::fixed<<std::setprecision(18)<<wham->prob[k][kj][kjl]<<std::endl;
            }//20240806
        }
    }

    fmin = 1.e200;
    for(i=0; i<wham->num_bins[0]; i++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                if(wham->prob[i][kj][kjl] > 1.e-200) {
                    wham->pmf[i][kj][kjl] = -1./beta*log(wham->prob[i][kj][kjl]);//20240806
                } else {
                    wham->pmf[i][kj][kjl] = 99999999.;
                }
                if(fmin > wham->pmf[i][kj][kjl]) fmin = wham->pmf[i][kj][kjl];
            }
        }
    }

    for(i=0; i<wham->num_bins[0]; i++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){//20240806
                wham->pmf[i][kj][kjl] -= fmin;
            }
        }
    }

    for(i=0; i<wham->num_bins[0]; i++) {
        for(kj=0; kj<wham->num_bins[1]; kj++) {
            for (kjl=0; kjl<wham->num_bins[2]; kjl++){
                x = (i + 0.5)*wham->bin_width[0] + wham->hist_min[0];
                y = (kj + 0.5)*wham->bin_width[1] + wham->hist_min[1];
                z = (kjl+ 0.5)*wham->bin_width[2] + wham->hist_min[2];//20240806
                //20240806
                std::cout<<"PMFF: "
                <<std::right<<std::setw(10)<<std::fixed<<std::setprecision(5)<<x
                <<std::right<<std::setw(10)<<std::fixed<<std::setprecision(5)<<y
                <<std::right<<std::setw(10)<<std::fixed<<std::setprecision(5)<<z
                <<std::right<<std::setw(32)<<std::fixed<<std::setprecision(18)<<wham->pmf[i][kj][kjl]
                <<std::right<<std::setw(32)<<std::fixed<<std::setprecision(18)<<wham->prob[i][kj][kjl]
                <<std::endl;;//20240806
            }
        }
    }
    
    std::cout<<"-----------------------------------------------------"<<std::endl;
    std::cout<<"bins: "<<std::endl;
    std::cout<<" "<<wham->hist_min[0]<<" "<<wham->hist_max[0]<<" "<<wham->num_bins[0]<<std::endl;
    std::cout<<" "<<wham->hist_min[1]<<" "<<wham->hist_max[1]<<" "<<wham->num_bins[1]<<std::endl;
    std::cout<<" "<<wham->hist_min[2]<<" "<<wham->hist_max[2]<<" "<<wham->num_bins[2]<<std::endl;
    std::cout<<"-----------------------------------------------------"<<std::endl;
    for (i=0;i<wham->num_windows;i++)
    {
        std::cout<<wham->F[i]<<std::endl;
    }
    
    for (i = 0;i<wham->num_bins[0];i++)
    {
        for(j = 0;j<wham->num_bins[1];j++)
        {
            delete histogram[i][j];
        }
        delete histogram[i];
    }
    delete histogram;
    std::cout<<"delete well"<<std::endl;

    for (i = 0 ;i<wham->num_windows;i++)
    {
        for (j=0;j<wham->nt[i];j++)
        {
            delete wham->rewX[i][j];
        }
        delete wham->rewX[i];
    }
    delete wham->rewX;
    std::cout<<"delete well"<<std::endl;

    for (i = 0 ;i<wham->num_windows;i++)
    {
        for (j=0;j<wham->nt[i];j++)
        {
            delete wham->X[i][j];
        }
        delete wham->X[i];
    }
    delete wham->X;
    std::cout<<"delete well"<<std::endl;

    for (i = 0 ;i<wham->num_windows;i++)
    {
        for (j=0;j<wham->nt[i];j++)
        {
            delete wham->ebw[i][j];
        }
        delete wham->ebw[i];
    }
    delete wham->ebw;
    std::cout<<"delete well"<<std::endl;

    for (i=0;i<wham->num_bins[0];i++)
    {        
        for (j=0;j<wham->num_bins[1];j++)
        {
            delete wham->prob[i][j];
            delete wham->pmf[i][j];
        }
        delete wham->prob[i];
        delete wham->pmf[i];
    }
    std::cout<<"delete well"<<std::endl;
    delete wham->prob;
    std::cout<<"delete well"<<std::endl;
    delete wham->pmf;
    std::cout<<"delete well"<<std::endl;
    delete [] wham->spring_constants;
    std::cout<<"delete well"<<std::endl;
    delete [] wham->bias_locations;
    std::cout<<"delete well"<<std::endl;
    // delete wham->dataFileName;
    // std::cout<<"delete well"<<std::endl;
    delete wham->nt;
    std::cout<<"delete well"<<std::endl;
    delete wham->partition;
    std::cout<<"delete well"<<std::endl;
    delete wham->kT;
    std::cout<<"delete well"<<std::endl;
    delete wham->F;
    std::cout<<"delete well"<<std::endl;
    delete wham;
    // std::cout<<"delete well"<<std::endl;
    // delete pmf;
    std::cout<<"delete well"<<std::endl;
    delete wham->ebf;
    std::cout<<"delete well"<<std::endl;

    return (1);
}

double calc_bias_pot(class wham_info *h, int index, double *coor, int BIAS)
{
    double springx,springy,springz;
    double dx,dy,dz;
    springx = h->spring_constants[index][0];
    springy = h->spring_constants[index][1];
    springz = h->spring_constants[index][2];
    dx = coor[0] - h->bias_locations[index][0];
    dy = coor[1] - h->bias_locations[index][1];
    dz = coor[2] - h->bias_locations[index][2];

    if (BIAS == 0) {
        dx = fabs(dx);
        if (dx > h->period[0]/2.)
        {
            dx -= h->period[0];
        }

        dy = fabs(dy);
        if (dy > h->period[1]/2.)
        {
            dy -= h->period[1];
        }

        dz = fabs(dz);//20240806
        if (dz > h->period[2]/2.)
        {
            dz -= h->period[2];
        }
        return 0.5*(springx*(dx*dx) + springy*(dy*dy) + springz*(dz*dz));
    } else if(BIAS == 1) {
        return ((springx*(1. - cos(dx)) + springy*(1. - cos(dy))));
    } else if(BIAS == 2) {
	    return (springx*dx);
    } else {
	    std::cout<<"unknown rstr type"<<BIAS<<" exit"<<std::endl;
        exit (1);
    }

    return (1);
}

int get_numwindows(char *file)
{
    std::string line;
    int num_windows;
    std::ifstream ifs;
    num_windows = 0;
    ifs.open(file,std::ios::in);
    while (getline(ifs,line))
    {
        if(!line.empty())
        {
            num_windows++;
        }
    }
    ifs.close();
    return num_windows;
}


#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
//a.out inputfile natom+2(nline) delta outfile
void comparen(int start,int end,int indexC,int &index,double **mat,double delta1,double delta2,int &acount);
void comparem(int start,int end,int indexC,int &index,double **mat,double mar[][3],double delta1,double delta2,int &acount,double &v1,double &v2);
void mostdeltaz(int start,int end,double **mat,double mar[][3],double deltaz,int &singal,int &indexRu);
// void comparen(int start,int end,int indexC,int &index,double **mat);
int main(int argc,char *argv[])
{
    double x,y,z,f_1,f_2,ene,ru_len,dist2d,tmpx,tmpy,tmpa,delta1,delta2;
    int a,n,natom,ncount,indexn,ind,acount,singal,indexRu,index,astep;
    double **fracmat,**fracmat33;
    double rmsd,rmsdz,rmsdsurf,rmsdsurfz,rmsdatom,rmsdzatom;
    double value1a,value2a,value3a,value4a,value5a,value6a;
    char elem[3],**element,**element33;
    int step,nline,nstep,ii,nv1,nv2;
    std::string line;
    std::ifstream ifs;
    std::ofstream ofs;

    if (argc < 7)
    {
        std::cout<<"./exec inputfile natom+2 delta1 delta2 index(temp) outfile\n";
        exit(0);
    }
    ifs.open(argv[1],std::ios::in);
    nline = atoi(argv[2]);
    delta1 = atof(argv[3]);
    delta2 = atof(argv[4]);
    index = atoi(argv[5]);
    ofs.open(argv[6],std::ios::out);
    n = nline;
    natom = nline - 2;
    // std::cout<<natom<<std::endl;
    ru_len = 2.7402961914621904*2.0;

    double matori[8][3] = {
        {1.58176 ,  1.06581 ,  6.41455},
        {3.95493  ,-0.30434   ,6.41455},
        {1.58176   ,3.80611   ,6.41455},
        {3.95493   ,2.43596  , 6.41455},
        {0.79117   ,2.43599  , 8.47202},
        {3.16433   ,1.06584  , 8.47202},
        {0.79117  ,-0.30431  , 8.47202},
        {3.16433  ,-1.67446  , 8.47202}
    };

    double matori33[13][3] = 
    {   {1.58176 ,  1.06581 ,  6.41455},//0
        {3.95493  ,-0.30434   ,6.41455},//1
        {1.58176   ,3.80611   ,6.41455},//2
        {3.95493   ,2.43596  , 6.41455},//3 //sub-ori
        {0.79117   ,2.43599  , 8.47202},//4 0 0.5 
        {3.16433   ,1.06584  , 8.47202},//5 0.5 0.5 
        {0.79117  ,-0.30431  , 8.47202},//6 0 0
        {3.16433  ,-1.67446  , 8.47202},//7 0.5 0//ori
        {0.79117  , 5.17628  , 8.47202},//8//top-new 0 1
        {5.53750  ,-3.04461  , 8.47202},//9 1 0
        {5.53750  , 2.43599  , 8.47202},//10 1 1 
        {3.16434  , 3.80613  , 8.47202},//11 0.5 0 -> 0.5 1
        {5.53750  ,-0.30431  , 8.47202}//12 0 0.5 -> 1 0.5  

    };

    element = new char *[natom];
    for (int i = 0; i<natom; i++)
    {
        element [i] = new char [3];
    }
    element33 = new char *[natom+5];
    for (int i = 0; i<natom+5; i++)
    {
        element33 [i] = new char [3];
    }
    fracmat = new double *[natom];
    for (int i = 0;i<natom;i++)
    {
        fracmat[i] = new double [3];
    }
    fracmat33 = new double *[natom+5];
    for (int i = 0;i<natom+5;i++)
    {
        fracmat33[i] = new double [3];
    }
    // std::cout<<"ALLOCATE ERROR!"<<std::endl;

    if (!ifs.is_open())
    {
        std::cout<<"No such file !"<<std::endl;
        exit(0);
    }
    nstep = 0;
    step = 0;
    value1a = 0.0;
    value2a = 0.0;
    value3a = 0.0;
    value4a = 0.0;
    value5a = 0.0;
    value6a = 0.0;

    ofs<<"step    RMSD    RMSD_Z     RMSD_L1    RMSD_L1Z   RMSDa    RMSDa_Z"<<std::endl;
    while (getline(ifs,line))
    {
        ncount = nstep % n;
        // std::cout<<ncount<<std::endl;
        std::stringstream linetemp(line);
        if (ncount == 0)
        {
            linetemp >> a;
            // ofs<<std::right<<std::setw(12)<<natom<<std::endl;
            // ofs<<std::right<<std::setw(12)<<natom+5<<std::endl;
        }else if(ncount == 1){
            linetemp >> ind >> ene;
            // ofs<<std::right<<std::setw(12)<<ind<<std::right<<std::setw(12)<<std::fixed<<std::setprecision(5)<<ene<<std::endl;
        }else {
            linetemp>>element[ncount-2]>>x>>y>>z;
            // element[ncount-2]=elem;
            // std::cout<<element[ncount-2]<<std::endl;
            fracmat[ncount-2][2]=z;
            x = x - 0.79117;
            y = y + 0.30431;
            f_1 = x/ru_len*2.0/sqrt(3.0);
            f_2 = y/ru_len+0.50*f_1;
            fracmat[ncount-2][0]=f_1;
            fracmat[ncount-2][1]=f_2;

            if (ncount == n-1){

                for (int i = 0;i<8;i++)//Ru surface at 2*2 area
                {
                    f_1 = fracmat[i][0];
                    f_2 = fracmat[i][1];
                    z = fracmat[i][2];
                    while (f_1 >= 1.0){//1.0 0.75
                        f_1 = f_1-1.0;
                    }
                    while (f_1 < -0.5){//-0.5 -0.25
                        f_1 = f_1+1.0;
                    }
                    while (f_2 >= 1.0){
                        f_2 = f_2-1.0;
                    }
                    while (f_2 < -0.5){
                        f_2 = f_2+1.0;
                    }
                    fracmat[i][0] = f_1;
                    fracmat[i][1] = f_2;
                    fracmat[i][2] = z;
                    element33[i] = element[i];
                    fracmat33[i][0] = fracmat[i][0];
                    fracmat33[i][1] = fracmat[i][1];
                    fracmat33[i][2] = fracmat[i][2];
                }

                for (int i = 8;i<13;i++)
                {
                    element33[i] = element[6];
                }

                fracmat33[8][0] = fracmat[6][0];
                fracmat33[8][1] = fracmat[6][1]+1.0;
                fracmat33[8][2] = fracmat[6][2];
                fracmat33[9][0] = fracmat[6][0]+1.0;
                fracmat33[9][1] = fracmat[6][1];
                fracmat33[9][2] = fracmat[6][2];
                fracmat33[10][0] = fracmat[6][0]+1.0;
                fracmat33[10][1] = fracmat[6][1]+1.0;
                fracmat33[10][2] = fracmat[6][2];
                //0.5 0 -> 0.5 1 extend 7
                fracmat33[11][0] = fracmat[7][0];
                fracmat33[11][1] = fracmat[7][1]+1.0;
                fracmat33[11][2] = fracmat[7][2];
                //0 0.5 -> 1 0.5 extend 4
                fracmat33[12][0] = fracmat[4][0]+1.0;
                fracmat33[12][1] = fracmat[4][1];
                fracmat33[12][2] = fracmat[4][2];
                // std::cout<<"well"<<std::endl;

                //co to 2*2 cell
                for (int i = 8;i<10;i++)
                {
                    f_1 = fracmat[i][0];
                    f_2 = fracmat[i][1];
                    z = fracmat[i][2];
                    while (f_1 >= 1.0){//1.0 0.75 0.5
                        f_1 = f_1-1.0;
                    }
                    while (f_1 < -0.5){//-0.5 -0.25 0.0
                        f_1 = f_1+1.0;
                    }
                    while (f_2 >= 1.0){
                        f_2 = f_2-1.0;
                    }
                    while (f_2 < -0.5){
                        f_2 = f_2+1.0;
                    }
                    fracmat[i][0] = f_1;
                    fracmat[i][1] = f_2;
                    fracmat[i][2] = z;
                    element33[i+5] = element[i];
                    fracmat33[i+5][0] = f_1;
                    fracmat33[i+5][1] = f_2;
                    fracmat33[i+5][2] = z;
                }

                // for (int i = 0;i<natom+5;i++)
                // {
                //     std::cout<<fracmat33[i][0]<<" "<<fracmat33[i][1]<<" "<<fracmat33[i][2]<<std::endl;
                // }
                // std::cout<<std::endl;

                double vx = 0.0;
                double dist = 0.0;
                for (int i = 0;i<2;i++)
                {
                    // vx += ((fracmat33[13][i]-fracmat33[14][i])*(fracmat33[13][i]-fracmat33[14][i]));
                    vx += pow((fracmat33[13][i]-fracmat33[14][i]),2);                     
                }
                // std::cout<<vx<<std::endl;
                dist = sqrt(vx);
                // std::cout<<dist<<std::endl;
                if (dist > 0.5)//0.5 0.25
                {
                    if ((fracmat33[13][0]-fracmat33[14][0]) > 0.50){fracmat33[14][0]=fracmat33[14][0]+1.0;}//1.0 0.5
                    if ((fracmat33[13][0]-fracmat33[14][0]) < -0.50){fracmat33[13][0]=fracmat33[13][0]+1.0;}
                    if ((fracmat33[13][1]-fracmat33[14][1]) > 0.50){fracmat33[14][1]=fracmat33[14][1]+1.0;}
                    if ((fracmat33[13][1]-fracmat33[14][1]) < -0.50){fracmat33[13][1]=fracmat33[13][1]+1.0;}
                }
                // for (int i = 0;i<natom;i++)
                // {
                //     std::cout<<std::fixed<<std::setprecision(5)<<fracmat[i][0]<<" "<<fracmat[i][1]<<" "<<fracmat[i][2]<<std::endl;
                // }
                // std::cout<<std::endl;
                //compare surface & co 

                for (int i = 0;i<natom+5;i++)
                {
                    f_1 = fracmat33[i][0];
                    f_2 = fracmat33[i][1];
                    x = f_1*ru_len*sqrt(3.0)/2.0;
                    y = (-0.50*f_1+f_2)*ru_len;
                    x = x + 0.79117;
                    y = y - 0.30431;
                    z = fracmat33[i][2];
                    fracmat33[i][0] = x;
                    fracmat33[i][1] = y;
                    fracmat33[i][2] = z;
                    // std::cout<<element[i]<<std::endl;
                    // ofs<<std::right<<std::setw(3)<<element33[i]
                    // <<std::right<<std::setw(14)<<std::fixed<<std::setprecision(5)<<x
                    // <<std::right<<std::setw(14)<<std::fixed<<std::setprecision(5)<<y
                    // <<std::right<<std::setw(14)<<std::fixed<<std::setprecision(5)<<z<<std::endl;
                    
                }

                // for (int i = 0;i<natom;i++)
                // {
                //     f_1 = fracmat[i][0];
                //     f_2 = fracmat[i][1];
                //     x = f_1*ru_len*sqrt(3.0)/2.0;
                //     y = (-0.50*f_1+f_2)*ru_len;
                //     x = x + 0.79117;
                //     y = y - 0.30431;
                //     z = fracmat[i][2];
                //     fracmat[i][0] = x;
                //     fracmat[i][1] = y;
                //     fracmat[i][2] = z;
                    
                // }


                // for (int i = 0;i<natom+5;i++)
                // {
                //     std::cout<<fracmat33[i][0]<<" "<<fracmat33[i][1]<<" "<<fracmat33[i][2]<<std::endl;
                // }
                // std::cout<<std::endl;

                //rmsd
                double value1 = 0.0;
                double value2 = 0.0;
                double value3 = 0.0;
                double value4 = 0.0;
                double value5 = 0.0;
                double value6 = 0.0;
                // std::cout<<value1<<std::endl;

                //rmsd all
                for (int k = 0;k<8;k++)
                {
                    for (int l = 0;l<3;l++)
                    {
                        // value1 += ((fracmat33[k][l]-matori33[k][l])*(fracmat33[k][l]-matori33[k][l]));
                        value1 += pow((fracmat33[k][l]-matori33[k][l]),2);
                    }
                    // value2 += ((fracmat33[k][2]-matori33[k][2])*(fracmat33[k][2]-matori33[k][2]));
                    value2 += pow((fracmat33[k][2]-matori33[k][2]),2);
                }
                // std::cout<<"value1 : "<<value1<<std::endl;
                //rmsd surf
                for (int k = 4;k<8;k++)
                {
                    for (int l = 0;l<3;l++)
                    {
                        // value3 += ((fracmat33[k][l]-matori33[k][l])*(fracmat33[k][l]-matori33[k][l]));
                        value3 += pow((fracmat33[k][l]-matori33[k][l]),2);
                    }
                    // value4 += ((fracmat33[k][2]-matori33[k][2])*(fracmat33[k][2]-matori33[k][2]));
                    value4 += pow((fracmat33[k][2]-matori33[k][2]),2);
                }
                // std::cout<<"value3 : "<<value3<<std::endl;
                // acount = 0;
                // comparen(4,13,13,indexn,fracmat33,delta1,delta2,acount);
                // // std::cout<<acount<<" "<<indexn<<std::endl;
                // // std::cout<<std::endl;
                // if (acount != 0)
                // {
                //     for (int l = 0; l < 3; l++)
                //     {
                //         // value5 +=  ((fracmat33[indexn][l]-matori33[indexn][l])*(fracmat33[indexn][l]-matori33[indexn][l]));
                //         value5 +=  pow((fracmat33[indexn][l]-matori33[indexn][l]),2);
                //     }
                //     // std::cout<<"value5 : "<<value5<<std::endl;
                //     // value6 = ((fracmat33[indexn][2]-matori33[indexn][2])*(fracmat33[indexn][2]-matori33[indexn][2]));
                //     value6 = pow((fracmat33[indexn][2]-matori33[indexn][2]),2);
                //     value5a += value5;
                //     value6a += value6;
                //     astep++;
                // }

                acount = 0;
                comparem(4,13,13,indexn,fracmat33,matori33,delta1,delta2,acount,value5,value6);
                if (acount != 0)
                {
                    value5a += sqrt(value5);
                    value6a += sqrt(value6);
                    astep++;
                }
                // std::cout<<indexn<<std::endl;
                // std::cout<<step<<std::endl;
                ofs<<std::right<<std::setw(10)<<step<<std::fixed<<std::setprecision(5)
                    <<std::right<<std::setw(10)<<sqrt(value1/8)
                    <<std::right<<std::setw(10)<<sqrt(value2/8)
                    <<std::right<<std::setw(10)<<sqrt(value3/4)
                    <<std::right<<std::setw(10)<<sqrt(value4/4)
                    <<std::right<<std::setw(10)<<sqrt(value5)
                    <<std::right<<std::setw(10)<<sqrt(value6)<<std::endl;
                value1a += sqrt(value1/8);
                value2a += sqrt(value2/8);
                value3a += sqrt(value3/4);
                value4a += sqrt(value4/4);           
                step++;               
            }
            
        }
        nstep++;
    }
    // std::cout<<value1a<<" "<<value2a<<" "<<value3a<<" "<<value4a<<" "<<value5a<<" "<<value6a<<std::endl;
    rmsd = value1a/step;
    rmsdz = (value2a/step);
    rmsdsurf = (value3a/step);
    rmsdsurfz = (value4a/step);
    rmsdatom = (value5a/astep);
    rmsdzatom = (value6a/astep);
    std::cout<<std::right<<std::setw(10)<<step<<std::fixed<<std::setprecision(5)
        <<std::right<<std::setw(10)<<rmsd
        <<std::right<<std::setw(10)<<rmsdz
        <<std::right<<std::setw(10)<<rmsdsurf
        <<std::right<<std::setw(10)<<rmsdsurfz
        <<std::right<<std::setw(10)<<rmsdatom
        <<std::right<<std::setw(10)<<rmsdzatom
        <<std::right<<std::setw(10)<<astep
        <<std::right<<std::setw(10)<<index<<std::endl;
    for (int i = 0 ;i < natom;i++)
    {
        delete fracmat[i];
    }
    delete fracmat;
    for (int i = 0 ;i < natom+5;i++)
    {
        delete fracmat33[i];
    }
    delete fracmat33;
    delete element;
    delete element33;
    ifs.close();
    // ofs.close();
    return 0;
}

// void comparen(int start,int end,int indexC,int &index,double **mat,double delta1,double delta2,int &acount)
// {
//     double *value,*r_value,tmp;
//     int inn = 0;
//     value = new double [end-start];
//     r_value = new double [end-start];
//     for (int i=0;i<end-start;i++)
//     {
//         tmp = 0.0;
//         for (int j=0;j<3;j++)
//         {
//             tmp += ((mat[i+start][j]-mat[indexC][j])*(mat[i+start][j]-mat[indexC][j]));
//         // std::cout<<value[i]<<std::endl;
//         }
//         value[i] = tmp;
//         r_value[i] = sqrt(value[i]);
//         // std::cout<<r_value[i]<<std::endl;
//         // std::cout<<std::endl;
//     }
//     for (int i=0;i<end-start;i++)
//     {            
//         if (r_value[i]>delta1 && r_value[i] < delta2)
//         {
//             acount++;
//             inn = i;
            
//         }
//         // std::cout<<r_value[i]<<std::endl;
//     }    
//     index = inn+start;
//     delete value;
//     delete r_value;
// }

void comparen(int start,int end,int indexC,int &index,double **mat,double delta1,double delta2,int &acount)
{
    double *value,*r_value,tmp,min;
    int inn = 0;
    int mininn = 0;
    value = new double [end-start];
    r_value = new double [end-start];
    for (int i=0;i<end-start;i++)
    {
        tmp = 0.0;
        for (int j=0;j<3;j++)
        {
            // tmp += ((mat[i+start][j]-mat[indexC][j])*(mat[i+start][j]-mat[indexC][j]));
            tmp += pow((mat[i+start][j]-mat[indexC][j]),2);
        // std::cout<<value[i]<<std::endl;
        }
        value[i] = tmp;
        r_value[i] = sqrt(value[i]);
        // std::cout<<r_value[i]<<std::endl;
        // std::cout<<std::endl;
    }
    min = 114514.1919810;
    for (int i = 0; i< end-start; i++)
    {
        if (r_value[i] < min)
        {
            min = r_value[i];
            mininn = i;
        }
    }
    // for (int i=0;i<end-start;i++)
    // {            
    //     if (r_value[i]>delta1 && r_value[i] < delta2)
    //     {
    //         acount++;
    //         inn = i;
            
    //     }
    //     std::cout<<r_value[i]<<std::endl;
    // }
    if(r_value[mininn] >= delta1 && r_value[mininn] < delta2)
    {
        acount++;
        // std::cout<<r_value[mininn]<<std::endl;
    }    
    
    index = mininn+start;
    // std::cout<<mininn<<" "<<index<<std::endl;
    delete value;
    delete r_value;
}

void comparem(int start,int end,int indexC,int &index,double **mat,double mar[][3],double delta1,double delta2,int &acount,double &v1,double &v2)
{
    double *value,*r_value,tmp,min;
    int inn = 0;
    int mininn = 0;
    value = new double [end-start];
    r_value = new double [end-start];
    for (int i=0;i<end-start;i++)
    {
        tmp = 0.0;
        for (int j=0;j<3;j++)
        {
            // tmp += ((mat[i+start][j]-mat[indexC][j])*(mat[i+start][j]-mat[indexC][j]));
            tmp += pow((mat[i+start][j]-mat[indexC][j]),2);
        // std::cout<<value[i]<<std::endl;
        }
        value[i] = tmp;
        r_value[i] = sqrt(value[i]);
        // std::cout<<r_value[i]<<std::endl;
        // std::cout<<std::endl;
    }
    // min = 114514.1919810;
    // for (int i = 0; i< end-start; i++)
    // {
    //     if (r_value[i] < min)
    //     {
    //         min = r_value[i];
    //         mininn = i;
    //     }
    // }
    for (int i=0;i<end-start;i++)
    {   
        double tmp1 = 0.0;
        double tmp2 = 0.0;         
        if (r_value[i] >= delta1 && r_value[i] < delta2)
        {
            acount++;
            inn = i;
            for (int j = 0; j < 3;j++)
            {
                tmp1 += pow((mat[i+start][j]-mar[i+start][j]),2);
            }
            tmp2 = pow((mat[i+start][2]-mar[i+start][2]),2);
            v1 += tmp1;
            v2 += tmp2;
            // std::cout<<tmp1<<" "<<tmp2<<std::endl;
        }
    //     std::cout<<r_value[i]<<std::endl;
    }
    v1 = v1/acount;
    v2 = v2/acount;

    index = mininn+start;
    // std::cout<<mininn<<" "<<index<<std::endl;
    delete value;
    delete r_value;
}

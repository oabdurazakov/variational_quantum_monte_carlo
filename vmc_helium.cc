#include<iostream>
#include<fstream>
#include<ctime>
#include<cstdlib>
#include<vector>
#include<cmath>
#include<algorithm>
#include<iomanip>

        using namespace std;

const int Nd = 3;
const int Nw = 500;
double a;
double b;
double dx = 0.2;
double energy;
double random_number(){return (double)rand()/(RAND_MAX + 1.0);}

double psi(double re1[], double re2[]){

	double r1 = 0, r2 = 0, r12 = 0;
	for (int i = 0; i < Nd; i++){
		r1 +=  re1[i] * re1[i];
		r2 +=  re2[i] * re2[i];
		r12 += (re1[i] - re2[i]) * (re1[i] - re2[i]);
	}
	r1 = sqrt(r1);
	r2 = sqrt(r2);
	r12 = sqrt(r12);
	double pow = -2.0 *(r1 + a * r1 * r1)/(1.0 + r1) - 2.0 * (r2 + a * r2 * r2)/(1.0 + r2) + 0.5 * r12 / (1.0 + b * r12); 
	return exp(pow);
}

double lap(double re1[], double re2[]) {

        double rp[Nd], rm[Nd], del2[Nd];
        double lap1 = 0.0, lap2 = 0.0;
        double h = 0.001;
        for (int i = 0; i < Nd; i++){
                for (int j = 0; j < Nd; j++){
                        if(j==i){
			rp[j] = re1[j] + h;
                        rm[j] = re1[j] - h;
			}
			else {
			rp[j] = re1[j];
                        rm[j] = re1[j];
			}
                }
                del2[i] = psi(rp, re2) - 2.0 * psi(re1,re2) + psi(rm,re2);
                lap1 += del2[i]/(h*h);
        }

        for (int i = 0; i < Nd; i++){
                for (int j = 0; j < Nd; j++){
                        if(j==i){
			rp[j] = re2[j] + h;
                        rm[j] = re2[j] - h;
			}
			else {
			rp[j] = re2[j];
                        rm[j] = re2[j];
			
			}
                }
                del2[i] = psi(re1, rp) - 2.0 * psi(re1,re2) + psi(re1,rm);
                lap2 += del2[i]/(h*h);
        }

        return (lap1 + lap2)/psi(re1,re2);
}



double eloc(double re1[], double re2[]){

	double r1 = 0, r2 = 0, r12 = 0;
	for (int i = 0; i < Nd; i++){
		r1 +=  re1[i] * re1[i];
		r2 +=  re2[i] * re2[i];
		r12 += (re1[i] - re2[i]) * (re1[i] - re2[i]);
	}
	r1 = sqrt(r1);
	r2 = sqrt(r2);
	r12 = sqrt(r12);


	double en = -0.5 * lap(re1,re2) - 2.0/r1 - 2.0/r2 + 1.0/r12;	
//	cout << "en:" << en << endl;
	return en;
}

struct Walker
{
	double r1[Nd], r2[Nd];
};

void initialize_walkers(vector<Walker> & walker){
	for (int i = 0; i < Nw; i++){
		for (int j = 0; j < Nd; j++){
			walker[i].r1[j] = random_number() - 0.5;
			walker[i].r2[j] = random_number() - 0.5;
		}
	}	
}

void metropolis_step(vector<Walker> & walker, int iw){
	double tempr1[Nd], tempr2[Nd];
	for (int j = 0; j < Nd; j++){
		tempr1[j]  = walker[iw].r1[j];
		tempr2[j]  = walker[iw].r2[j];
		tempr1[j] += dx * (2.0 * random_number() - 1.0);
		tempr2[j] += dx * (2.0 * random_number() - 1.0);
	}
	double prp = psi(tempr1,tempr2) * psi(tempr1,tempr2);
	double pr  = psi(walker[iw].r1, walker[iw].r2) * psi(walker[iw].r1, walker[iw].r2);
	if(prp/pr >= random_number()){
		for (int j = 0; j < Nd; j++){
			walker[iw].r1[j] = tempr1[j]; 
			walker[iw].r2[j] = tempr2[j];
		}
	}
	energy += eloc(walker[iw].r1, walker[iw].r2);	
} 

main(){
	srand(time(NULL));
	vector<Walker> walker(Nw);
	initialize_walkers(walker);
	double Nther = 2000;
	ofstream outfile("energy.dat");
	a = 0.0;	
	while(a<=1.5){
		b = 0.0;	
		while(b<=1.5){	
			energy = 0.0;
			for (int j = 0; j < Nther; j++){
				for (int i = 0; i < Nw; i++){metropolis_step(walker,i);}
			}
			outfile << a << " "  << b << " " << energy/Nw/Nther << endl;
			b += 0.05;
	}
		a += 0.05;
	}
	
}

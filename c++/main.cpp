#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include "eigen-3.4-rc1/Eigen/Dense"
#include <fstream>
#include <cstdio>

using namespace Eigen;
using namespace std;
 

MatrixXf lk();
VectorXf FE(float nelx, float nely,MatrixXf x, float penal);
MatrixXf filter(float nelx,float nely,float rmin, MatrixXf x, MatrixXf dc);
MatrixXf OC(float nelx,float nely, MatrixXf x,float volfrac, MatrixXf dc);
void top(float nelx,float nely, float volfrac, float penal, float rmin);


int main(){
    
    top(30,10,0.5,3.0,1.5);
      
    return 0;
}

void top(float nelx,float nely, float volfrac, float penal, float rmin){
     // INITIALIZE
    MatrixXf x = MatrixXf::Ones(nely,nelx)*volfrac; // density matrix
    int loop = 0;
    float change = 1;
    MatrixXf dc = MatrixXf::Ones(nely,nelx); // Sensitivity matrix
    //start iteration
    while (change > 0.01){
        loop = loop+1;
        MatrixXf xold = x;
    // FE-ANALYSIS
        VectorXf U = FE(nelx,nely,x,penal);
        // cout<< U;
    // OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS 
        MatrixXf KE = lk();
        float c = 0;
        for (int ely = 0; ely < nely; ely++){
            for (int elx = 0; elx < nelx; elx++){
                int n1 = (nely+1)*(elx)+ely+1;
                int n2 = (nely+1)* (elx+1)   +ely+1;
                vector <int> ind  {2*n1-2,2*n1-1, 2*n2-2,2*n2-1, 2*n2,2*n2+1, 2*n1,2*n1+1};
                VectorXf Ue = U(ind);
                c = c + pow(x(ely,elx),penal)*float((Ue.transpose()*KE)*Ue);
                dc(ely,elx) = -penal*(pow(x(ely,elx),(penal-1)))*float((Ue.transpose()*KE)*Ue);
            }
        }
        // FILTERING THE SENSITIVITIES
        // cout<< dc <<endl;
        dc = filter(nelx,nely,rmin,x,dc);
        // UPDATE DESIGN BY THE OPTIMALITY CRITERIA METHOD
        x = OC(nelx,nely,x,volfrac,dc);      
        // PRINT RESULTS
        change = ((x.array()-xold.array()).abs()).maxCoeff();

        cout<<"It.: "<< loop << "  Obj.: " << c << "  vol:. "<< x.sum()/(nelx*nely)<<  "  ch.: " <<change<< endl;

            char buffer[15];

        sprintf(buffer, "image%d.ppm", loop);
        ofstream image;
        image.open(buffer);
        if(image.is_open()){
            image << "P1"<< endl;
            image << nelx << "  "<< nely << endl;
            // image << "255" << endl;

    
            for (int ely = 0; ely < nely; ely++){
                for(int elx = 0; elx < nelx; elx++){
                    image << round(x(ely,elx))<< " ";
                }
                image<<endl;
            }
        }

    }

    cout<<x;
}

////////////////// ELEMENT STIFFNESS MATRIX ///////////////////
MatrixXf lk(){
    float E = 1; 
    float nu = 0.3;
    VectorXf k(9);
    k << 0, 1.0/2-nu/6,1.0/8+nu/8,-1.0/4-nu/12,-1.0/8+3*nu/8,-1.0/4+nu/12,-1.0/8-nu/8,nu/6,1.0/8-3*nu/8;
    MatrixXf r(8,8);
    r << k(1), k(2), k(3), k(4), k(5), k(6), k(7), k(8),
         k(2), k(1), k(8), k(7), k(6), k(5), k(4), k(3),
         k(3), k(8), k(1), k(6), k(7), k(4), k(5), k(2),
         k(4), k(7), k(6), k(1), k(8), k(3), k(2), k(5),
         k(5), k(6), k(7), k(8), k(1), k(2), k(3), k(4),
         k(6), k(5), k(4), k(3), k(2), k(1), k(8), k(7),
         k(7), k(4), k(5), k(2), k(3), k(8), k(1), k(6),
         k(8), k(3), k(2), k(5), k(4), k(7), k(6), k(1);
    MatrixXf KE = (E/(1-pow(nu,2)))*r; //  Element Stiffness matrix
    
    return KE;
}
//////////////////////// FE-ANALYSIS //////////////////////////////////
VectorXf FE(float nelx, float nely,MatrixXf x, float penal){
    MatrixXf KE = lk();
    MatrixXf K = MatrixXf::Zero(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
    VectorXf F = VectorXf::Constant(2*(nelx+1)*(nely+1),0);
    VectorXf U = VectorXf::Constant(2*(nelx+1)*(nely+1),0);
    for (int elx =0; elx < nelx; elx++){
        for (int ely = 0; ely < nely; ely++){
            float n1 = ((nely+1)*(elx))+ely+1;
            float n2 = (nely+1)* (elx+1)   +ely+1;
            vector <float> edof = {2*n1-2,2*n1-1, 2*n2-2,2*n2-1, 2*n2,2*n2+1, 2*n1,2*n1+1};
            K(edof,edof) = K(edof,edof).array() + (pow(x(ely,elx),penal)*KE).array();
            
        }
    }
    // DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
    
    F(1,0) = -1;
    vector <int> alldofs;
    for (int i = 0; i<(2*(nely+1)*(nelx+1));i++){alldofs.push_back(i);}
    vector <int> fixeddofs;
    for (int i = 0; i<(2*(nely+1)) ;i++){
        if (i%2==0) {fixeddofs.push_back(i);}
    }
    fixeddofs.push_back(2*(nelx+1)*(nely+1)-1);
    
    sort(alldofs.begin(), alldofs.end());
    sort(fixeddofs.begin(), fixeddofs.end());
    vector <int> freedofs;
    set_symmetric_difference(alldofs.begin(), alldofs.end(),fixeddofs.begin(), fixeddofs.end(),back_inserter(freedofs));
    // SOLVING
    U(freedofs) = (K(freedofs,freedofs).inverse())*F(freedofs);
    U(fixeddofs)= VectorXf::Constant(fixeddofs.size(),0);
    return U;
}
//////////////////// MESH-INDEPENDENCY FILTER ///////////////////////
MatrixXf filter(float nelx,float nely,float rmin, MatrixXf x, MatrixXf dc){
    
    MatrixXf dcn = MatrixXf::Zero(nely,nelx);
    for (int i = 0; i < nelx; i++){
        for (int j = 0; j < nely; j++){
            float sum = 0.0;
            for (int k= max((i+1)-floor(rmin),float(1)); k < min((i+1)+floor(rmin),nelx)+1; k++){
                for (int l =max((j+1)-floor(rmin),float(1)); l < min((j+1)+floor(rmin), nely)+1; l++){
                float fac = rmin-sqrt(pow((i+1)-k,2) + pow((j+1)-l,2));
                sum = sum+max(float(0),fac);
                dcn(j,i) = dcn(j,i)+ max(float(0),fac)*x(l-1,k-1)* dc(l-1,k-1);
                }
            }
            dcn(j,i) = dcn(j,i) / (x(j,i)*sum);
        }
    }
    return dcn;
}
/////////////////// OPTIMALITY CRITERIA UPDATE /////////////
MatrixXf OC(float nelx,float nely, MatrixXf x,float volfrac, MatrixXf dc){
    float l1 = 0;
    float l2 = 100000;
    float move = 0.2;
    float lmid;
    MatrixXf xnew = MatrixXf::Zero(nely,nelx);

    while ((l2-l1) > 0.0001){
        
        lmid = 0.5*(l2+l1);
        MatrixXf temp1 = MatrixXf::Ones(nely,nelx);
        MatrixXf temp2 = MatrixXf::Ones(nely,nelx)*0.001;
        xnew = (temp2.array().max((x.array()-move).max(temp1.array().min((x.array()+move).min(x.array()*(sqrt((-dc/lmid).array())))))));
        if ((xnew.sum() - volfrac*nelx*nely) >= 0){
            l1 = lmid;
        }else{
            l2 = lmid;
        }
    }
    return xnew;
}
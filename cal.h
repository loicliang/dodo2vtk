#ifndef H_cal
#define H_cal


#include <math.h>
#include <vector>
#include <limits>
#define alphaDV 0.403226
#define stabilized_cycle 5
const double SIGMA=152.5;
struct tensor1;
struct tensor2;
double maxdouble=(std::numeric_limits<double>::max)(); 
double mindouble=(std::numeric_limits<double>::min)();


double sigmah(double t)
{
    double temp;
    t=t-floor(t);
    if(t>=0.0 && t<=0.25) temp=4.0*t*SIGMA;
        else if(t>0.25 && t<=0.75) temp=-4.0*SIGMA*(t-0.5);
        else if(t>0.75 && t<1.0)  temp=4.0*SIGMA*(t-1.0);
        else temp=0.0;
    return temp;
};
struct array12
{
    double a[12];
    void operator = (const double t[12])
    {
        for (int i=0;i<12;i++)
        {
            this->a[i]=t[i];
        }
    }
};


struct Elt
{
  int gn=0;
  int num=0;
  bool f=false;
  double volume=0.0;
};
struct Grain
{
    int num=0;
    int elts=0;
    bool f=false;
    double volume=0.0;
    std::vector<int> n_elt;
};



 struct tensor1
{
    double a1,a2,a3;

    double operator * (const tensor1 &A)
    {   double temp;
        temp=this->a1*A.a1+this->a2*A.a2+this->a3*A.a3;
        return temp;
    }
    tensor1 operator * (const double &A)
    {   tensor1 temp;
        temp.a1=this->a1*A;
        temp.a2=this->a2*A;
        temp.a3=this->a3*A;
        return temp;
    }
    tensor1 operator + (const tensor1 &A)
    {   tensor1 temp;
        temp.a1=this->a1+A.a1;
        temp.a2=this->a2+A.a2;
        temp.a3=this->a3+A.a3;
        return temp;
    }
    tensor1 operator - (const tensor1 &A)
    {   tensor1 temp;
        temp.a1=this->a1-A.a1;
        temp.a2=this->a2-A.a2;
        temp.a3=this->a3-A.a3;
        return temp;
    }


    void operator = (const double A[3])
    {

        this->a1=A[0];
        this->a2=A[1];
        this->a3=A[2];
    }

    double mode (void)
    {   double temp;
        temp=sqrt(this->a1*this->a1+this->a2*this->a2+this->a3*this->a3);
        return temp;
    }

};

struct tensor2
{
    double a11,a12,a13,a21,a22,a23,a31,a32,a33;

    tensor1 operator * (const tensor1 &A)
    {   tensor1 temp;
        temp.a1=this->a11*A.a1+this->a12*A.a2+this->a13*A.a3;
        temp.a2=this->a21*A.a1+this->a22*A.a2+this->a23*A.a3;
        temp.a3=this->a31*A.a1+this->a32*A.a2+this->a33*A.a3;
        return temp;
    }
    tensor2 operator + (const tensor2 &A)
    {
        tensor2 temp;
        temp.a11=this->a11+A.a11;
        temp.a12=this->a12+A.a12;
        temp.a13=this->a13+A.a13;
        temp.a21=this->a21+A.a21;
        temp.a22=this->a22+A.a22;
        temp.a23=this->a23+A.a23;
        temp.a31=this->a31+A.a31;
        temp.a32=this->a32+A.a32;
        temp.a33=this->a33+A.a33;
        return temp;
    }
    tensor2 operator / (const double &A)
    {
        tensor2 temp;
        temp.a11=this->a11/A;
        temp.a12=this->a12/A;
        temp.a13=this->a13/A;
        temp.a21=this->a21/A;
        temp.a22=this->a22/A;
        temp.a23=this->a23/A;
        temp.a31=this->a31/A;
        temp.a32=this->a32/A;
        temp.a33=this->a33/A;
        return temp;
    }
    tensor2 operator * (const double &A)
    {
        tensor2 temp;
        temp.a11=this->a11*A;
        temp.a12=this->a12*A;
        temp.a13=this->a13*A;
        temp.a21=this->a21*A;
        temp.a22=this->a22*A;
        temp.a23=this->a23*A;
        temp.a31=this->a31*A;
        temp.a32=this->a32*A;
        temp.a33=this->a33*A;
        return temp;
    }    
        
    void operator = (const double A[6])
    {

        this->a11=A[0];
        this->a22=A[1];
        this->a33=A[2];
        this->a21=A[3];
        this->a12=A[3];
        this->a23=A[4];
        this->a32=A[4];
        this->a13=A[5];
        this->a31=A[5];

    }

};



 inline   tensor1  vectorproduct(const tensor1 &A,const tensor1 &B)
    {
        tensor1 temp;
        temp.a1=A.a2*B.a3-A.a3*B.a2;
        temp.a2=A.a3*B.a1-A.a1*B.a3;
        temp.a3=A.a1*B.a2-A.a2*B.a1;
        return temp;
    };

 inline   tensor2 kroneckerproduct(const tensor1 &A, const tensor1 &B)
    {
        tensor2 temp;
        temp.a11=A.a1*B.a1;
        temp.a12=A.a1*B.a2;
        temp.a13=A.a1*B.a3;
        temp.a21=A.a2*B.a1;
        temp.a22=A.a2*B.a2;
        temp.a23=A.a2*B.a3;
        temp.a31=A.a3*B.a1;
        temp.a32=A.a3*B.a2;
        temp.a33=A.a3*B.a3;
        return temp;
    };

inline    double doubleinnerproduct(const tensor2 &A, const tensor2 &B)
    {
      double temp;
      temp=A.a11*B.a11+A.a12*B.a21+A.a13*B.a31+A.a21*B.a12+A.a22*B.a22+A.a23*B.a32+A.a31*B.a13+A.a32*B.a23+A.a33*B.a33;
      return temp;
    };


 inline   double vonmises(const tensor2 &A)
    {
        double temp;
        temp= sqrt( ( (A.a11-A.a22)*(A.a11-A.a22) + (A.a22-A.a33)*(A.a22-A.a33) + (A.a33-A.a11)*(A.a33-A.a11) + 6.0*( A.a12*A.a12 + A.a23*A.a23+A.a31*A.a31) )*0.5) ;
        return temp;
    };



#endif

#include <stdio.h>
#include <math.h>

void Taylor(int b, long n, double epsilon, double an, double sum, double sum1, int i, double* xvalues, double* yvalues, int col, int kol, double* xvaluesCheb, double* yvaluesCheb);
double InterpolateLagrangePolynomial (double x, double* x_values, double* y_values, int size);
void Otrisovka (double l, int mashtab);
double Fx (double epsilon, double x); 
double Newton(double x, int n, double* x_arr, double* y_arr);

int main() {
    double an=0, sum=0, sum1=0;
    int i=0;
    int col=10;
    long n=0;
    double epsilon = pow(10, -4);
    double a=0, b=2;
    double yvalues[100];
    double xvalues[100];
    double xvaluesCheb[100];
    double yvaluesCheb[100];
    printf("Task1.\n");
    for (int kol=6; kol<90; kol++) {
        //printf("%d\n", kol);
        Taylor(b, n, epsilon, an, sum, sum1, i, xvalues, yvalues, col, kol, xvaluesCheb, yvaluesCheb);
    }
    printf("\n\n");
    printf("Task2.\n");
    return 0;
}
    
double Newton(double x, int n, double* x_arr, double* y_arr){
    
    double sum = y_arr[0];
    for(int i = 1; i < n; ++i){
        
        double F = 0;
        for(int j = 0; j <= i; ++j){
            
            double den = 1;
            for(int k = 0; k <= i; ++k)
                if (k != j) 
                    den *= (x_arr[j] - x_arr[k]);
            F += y_arr[j]/den;
        }
        
        for(int k = 0; k < i; ++k)
            F *= (x - x_arr[k]);
        sum += F;
    }
    return sum;
}

double Fx (double epsilon, double x) {
    double sum = 0;
    int n = 0;
    double an = x;
    double q=0;
        while (fabs(an) > epsilon) {
            sum += an;
            q = -(x * x * (2 * n + 1)) / ((n + 1) * (2 * n + 3));
            an *= q;
            n++;
        }
        sum *= 2 / sqrt(M_PI);
    return sum;
}

void Taylor (int b, long n, double epsilon, double an, double sum, double sum1, int i, double* xvalues, double* yvalues, int col, int kol, double* xvaluesCheb, double* yvaluesCheb) {
    int a=0;
    double max=-100;
    double h = (double)(b - a) / (kol-1);
    double h1 = (double)(b - a) / 10;

    for (int g = 0; g <kol; g += 1) {
        double x=(b + a) / 2 + (b - a) / 2 * cos((2 * g + 1) * M_PI /(2 * kol + 2));
        xvaluesCheb[g]=x;
        //printf("%lf %d\n", xvaluesCheb[g], kol);
    }

    for (int g = 0; g <=kol; g += 1) {
        double x=(b + a) / 2 + (b - a) / 2 * cos((2 * g + 1) * M_PI /(2 * kol + 2));
        yvaluesCheb[g]=Fx(epsilon, x);
        //printf("%.1lf %lf\n", xvaluesCheb[g], yvaluesCheb[g]);
    }

    for (int top=0; top<kol; top++){
        double xi=a+top*h;
        xvalues[top]=xi;
        //printf("%lf\n", xvalues[top]);
        }
    //printf("%d ", kol);
    //printf(" x    f(x)\n");
    for (i = 0; i <kol; i += 1) {
        double x=a+i*h;
        yvalues[i]=Fx(epsilon, x);
        //printf("%lf\n", yvalues[i]);
    }  

    for (i = 0; i <=10; i += 1) {
        double x=a+i*h1;
        double xCheb=(b + a) / 2 + (b - a) / 2 * cos((2 * i + 1) * M_PI /(2 * 10 + 2));
        double Lagrange=InterpolateLagrangePolynomial(x, xvalues, yvalues, kol);
        double LangrangeCheb=InterpolateLagrangePolynomial(xCheb, xvaluesCheb, yvaluesCheb, kol);
        double newtown=Newton(x, kol, xvaluesCheb, yvaluesCheb);
        //printf ("%lf\n", newtown);
        //printf("%.1lf %.9lf\n", xCheb, fabs(LangrangeCheb));

        if (max<fabs(LangrangeCheb-Fx(epsilon, xCheb))) {
            max=fabs(LangrangeCheb-Fx(epsilon, xCheb));
        }

        //if (max<fabs(Lagrange-Fx(epsilon, x))) {
          //  max=fabs(Lagrange-Fx(epsilon, x));
        //}
    }

    printf("%.9lf\n", max);
    return;
}

double InterpolateLagrangePolynomial (double x, double* x_values, double* y_values, int size)
{
	double lagrange_pol = 0;
	double basics_pol;

	for (int i = 0; i < size; i++)
	{
		basics_pol = 1;
		for (int j = 0; j < size; j++)
		{
			if (j == i) continue;
			basics_pol *= (x - x_values[j])/(x_values[i] - x_values[j]);		
		}
		lagrange_pol += basics_pol*y_values[i];
	}
	return lagrange_pol;
}


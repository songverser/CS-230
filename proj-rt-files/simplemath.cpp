#include "simplemath.h"

sol QuadEq(double a, double b, double c)
{
    float x1, x2, discriminant;//, realPart, imaginaryPart;
    discriminant = b * b - 4 * a * c;
    struct sol ans;
    if (discriminant >= 0)
    {
        x1 = (-b + sqrt(discriminant)) / (2 * a);
        x2 = (-b - sqrt(discriminant)) / (2 * a);
        // std::cout << "Roots are real." << std::endl;
        // std::cout << "x1 = " << x1 << std::endl;
        // std::cout << "x2 = " << x2 << std::endl;
        ans.solvable = true;
        ans.solution[0] = x1;
        ans.solution[1] = x2;
    }
    else
    {
        // realPart = -b / (2 * a);
        // imaginaryPart = sqrt(-discriminant) / (2 * a);
        // std::cout << "No real solution" << std::endl;
        // std::cout << "x1 = " << realPart << "+" << imaginaryPart << "i" << std::endl;
        // std::cout << "x2 = " << realPart << "-" << imaginaryPart << "i" << std::endl;
        ans.solvable = false;
    }

    return ans;
}

double minpositive(double a, double b){
    if (a<b && a>small_t) return a;
    else if(b>small_t) return b;
    else return 0;
}

int gettreelevel(int i){
    int ans = 0;
    while(pow(2,ans)-2<i){
        ans++;
    }
    return ans;
}

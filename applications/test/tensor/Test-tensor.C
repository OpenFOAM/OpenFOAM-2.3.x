#include "tensor.H"
#include "symmTensor.H"
#include "transform.H"
#include "stringList.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    tensor t1(1, 2, 3, 4, 5, 6, 7, 8, 9);
    tensor t2(1, 2, 3, 1, 2, 3, 1, 2, 3);

    tensor t3 = t1 + t2;

    Info<< t3 << endl;

    tensor t4(3,-2,1,-2,2,0,1, 0, 4);

    Info<< inv(t4) << endl;
    Info<< (inv(t4) & t4) << endl;

    Info<< t1.x() << t1.y() << t1.z() << endl;

    tensor t6(1,0,-4,0,5,4,-4,4,3);
    //tensor t6(1,2,0,2,5,0,0,0,0);
    Info<< "tensor " << t6 << endl;
    vector e = eigenValues(t6);
    Info<< "eigenvalues " << e << endl;
    tensor ev = eigenVectors(t6);
    Info<< "eigenvectors " << ev << endl;

    Info<< "Check determinant " << e.x()*e.y()*e.z() << " " << det(t6) << endl;

    Info<< "Check eigenvectors "
        << (eigenVector(t6, e[0]) & t6) << e[0]*eigenVector(t6, e[0]) << " "
        << (eigenVector(t6, e[1]) & t6) << e[1]*eigenVector(t6, e[1]) << " "
        << (eigenVector(t6, e[2]) & t6) << e[2]*eigenVector(t6, e[2])
        << endl;

    Info<< "Check eigenvalues for symmTensor "
        << eigenValues(symm(t6)) - eigenValues(tensor(symm(t6))) << endl;

    Info<< "Check eigenvectors for symmTensor "
        << eigenVectors(symm(t6)) - eigenVectors(tensor(symm(t6))) << endl;

    tensor t7(1, 2, 3, 2, 4, 5, 3, 5, 6);

    Info<< "Check transformation "
        << (t1 & t7 & t1.T()) << " " << transform(t1, t7) << endl;

    symmTensor st1(1, 2, 3, 4, 5, 6);
    symmTensor st2(7, 8, 9, 10, 11, 12);

    Info<< "Check symmetric transformation "
        << transform(t1, st1) << endl;

    Info<< "Check for dot product of symmetric tensors "
        << (st1 & st2) << endl;

    vector v1(1, 2, 3);

    Info<< sqr(v1) << endl;
    Info<< symm(t7) << endl;
    Info<< twoSymm(t7) << endl;
    Info<< magSqr(st1) << endl;
    Info<< mag(st1) << endl;

    Info<< (symm(t7) && t7) - (0.5*(t7 + t7.T()) && t7) << endl;
    Info<< (t7 && symm(t7)) - (t7 && 0.5*(t7 + t7.T())) << endl;

    return 0;
}

#include <iostream>
#include "mpreal.h"
#include <vector>
#include <random>
#include "CMA.hpp"
#include <cmaes.h>

using namespace std;
using namespace mpfr;
using namespace libcmaes;

class HyperTrig {
public:
    vector<mpreal> a, b;
    mpreal c;
    mpreal omega;

    HyperTrig(int n) {
        a.resize(n);
        b.resize(n);
    }
    void fromCosSquare (mpreal om, const vector<mpreal>& phi);
    mpreal cv2();
};

static mpreal mpOne (1);
static int precision;

void compCoeffs(const vector<mpreal>& phi,
                HyperTrig& target,
                HyperTrig& buffer) {

    target.a.resize(phi.size());
    target.b.resize(phi.size());
    for (int i=0; i<phi.size(); i++) {
        if (i==0) {
            target.a[0] = cos(phi[0])/2;
            target.b[0] = sin(phi[0])/2;
            target.c = mpOne/2;
        }
        else if (i==1) {
            target.a[0] = (cos(phi[0])+cos(phi[1])) / 4;
            target.a[1] = cos(phi[0]+phi[1]) / 8;
            target.b[0] = (sin(phi[0])+sin(phi[1])) / 4;
            target.b[1] = sin(phi[0]+phi[1]) / 8;
            target.c = cos(phi[0]-phi[1])/8 + mpOne/4;
        }
        else {
            buffer = target;
            mpreal cosN = cos(phi[i]);
            mpreal sinN = sin(phi[i]);
            target.a[0] = buffer.a[0]/2 + buffer.a[1]*cosN/4 + buffer.b[1]*sinN/4 + buffer.c*cosN/2;
            target.b[0] = buffer.b[0]/2 - buffer.a[1]*sinN/4 + buffer.b[1]*cosN/4 + buffer.c*sinN/2;
            for (int k=1; k<i-1; k++) {
                target.a[k] = buffer.a[k-1]*cosN/4 - buffer.b[k-1]*sinN/4 + buffer.a[k]/2 + buffer.a[k+1]*cosN/4 + buffer.b[k+1]*sinN/4;
                target.b[k] = buffer.b[k-1]*cosN/4 + buffer.a[k-1]*sinN/4 + buffer.b[k]/2 - buffer.a[k+1]*sinN/4 + buffer.b[k+1]*cosN/4;
            }
            target.a[i-1] = buffer.a[i-2]*cosN/4 - buffer.b[i-2]*sinN/4 + buffer.a[i-1]/2;
            target.b[i-1] = buffer.b[i-2]*cosN/4 + buffer.a[i-2]*sinN/4 + buffer.b[i-1]/2;
            target.a[i-0] = buffer.a[i-1]*cosN/4 - buffer.b[i-1]*sinN/4;
            target.b[i-0] = buffer.b[i-1]*cosN/4 + buffer.a[i-1]*sinN/4;
            target.c = buffer.c/2 + buffer.a[0]*cosN/4 + buffer.b[0]*sinN/4;
        }
    }
}

void HyperTrig::fromCosSquare(mpreal om, const vector<mpreal>& phi) {
    omega = om;
    HyperTrig buffer(200);
    compCoeffs(phi, *this, buffer);
}

mpreal HyperTrig::cv2() {

    mpreal m0 = c;
    mpreal m1 = -c;
    mpreal m2 = c*2;
    for (int i=0; i<a.size(); i++) {
        mpreal x =  mpOne + omega*omega*(i+1)*(i+1);
        m0 += (a[i] + b[i]*(i+1)*omega) / x;
        m1 += (a[i]*(x-mpreal(2)) - mpreal(2)*b[i]*(i+1)*omega) / x / x;
        m2 += (a[i]*(mpreal(8)-x*6) - mpreal(2)*b[i]*(i+1)*omega*(x-mpreal(4))) / x / x / x;
    }
    return m2*m0/m1/m1-mpOne;
}

pair<mpreal, vector<mpreal>> rechenberg_optimize(mpreal initialOmega, const vector<mpreal>& initialPhi, int trials=5) {

    int N = initialPhi.size();
    HyperTrig target(N);

    vector<mpreal> optPhi = initialPhi;
    mpreal optOmega = initialOmega;
    target.fromCosSquare(optOmega, optPhi);
    mpreal optCv2 = target.cv2();

    vector<mpreal> bestPhi = optPhi;
    mpreal bestOmega = optOmega;
    mpreal bestCv2 = optCv2;

    vector<mpreal> yPhi (N);

    random_device rd;
    mt19937 gen (rd());
    std::normal_distribution<> d (0.0, 1.0);

    for (int k=0; k<trials; k++) {

        optPhi = initialPhi;
        optOmega = initialOmega;
        target.fromCosSquare(optOmega, optPhi);
        optCv2 = target.cv2();

        mpreal sigma (1);
        mpreal c (0.95);
        int success = 0;

        for (int i=0; sigma>1e-6; i++) {

            mpreal yOmega = abs(optOmega + sigma*d(gen));
//            yOmega = initialOmega;
            for (int j=0; j<N; j++)
                yPhi[j] = optPhi[j] + sigma*d(gen);

            target.fromCosSquare(yOmega, yPhi);

            mpreal yCv2 = target.cv2();
            if (yCv2 < optCv2) {
                success++;
                optCv2 = yCv2;
                optOmega = yOmega;
                optPhi = yPhi;
            }

            if (i%20 == 0) {
                if (success < 5)
                    sigma *= c;
                else
                    sigma /= c;
                success = 0;
            }
            if (i%5000==0 && (i-1)%5000!=0) {
                cout << "* " << i << ", " << i << ", cv2: " << optCv2 << ", sigma: " << sigma << endl;
                cout << "* omega = " << optOmega << endl;
                cout << "* phi = {";
                for (int i=0; i<N; i++)
                    if (i<N-1)
                        cout << optPhi[i] << ",";
                    else
                        cout << optPhi[i] << "}" << endl;
            }
        }
        cout << "** trial " << k+1 << ": cv2(" << N << ") = " << optCv2 << endl;
        if (optCv2 < bestCv2) {
            bestCv2 = optCv2;
            bestOmega = optOmega;
            bestPhi = optPhi;
        }
    }
    cout << "cv2N" << N << "=" << bestCv2 << endl;
    cout << "omegaN" << N << "=" << bestOmega << endl;
    cout << "phiN" << N << "={";
    for (int i=0; i<N; i++)
        if (i<N-1)
            cout << bestPhi[i] << ",";
        else
            cout << bestPhi[i] << "}" << endl;

//    target.fromCosSquare(optOmega, optPhi);
    return make_pair (bestOmega, bestPhi);
}

mpreal cv2 (vector<mpreal> u) {

    HyperTrig target(u.size()-1);

    vector<mpreal> phi (u.begin()+1, u.end());
    mpreal omega = mpfr::abs(u[0]);
    target.fromCosSquare(omega, phi);
    mpreal cv2 = target.cv2();
    return cv2;
}

double funWithGap(double omega, double place, double gap, int N) {
    double delta = (2.0*M_PI-gap) / N;
    int i = floor(place/delta+0.5);
    vector<mpreal> params;
    params.push_back(omega);
    for (int j=1; j<=i; j++)
        params.push_back((j-0.5)*delta-M_PI);
    for (int j=i+1; j<=N; j++)
        params.push_back((j-0.5)*delta+gap-M_PI);
    return cv2(params).toDouble();
}

pair<mpreal,vector<mpreal>> initialGuess (int N) {

    double sigma = 1.0;
    double c = 0.9;
    int success = 0;

    double place = 4.0;
    double gap = 2.0;
    double omega = 0.5;

    double bestVal = 1e12;
    random_device rd;
    mt19937 gen (rd());
    std::normal_distribution<> d (0.0, 1.0);

    for (int i=0; sigma>1e-5; i++) {

        double thisPlace = place + sigma*d(gen);
        double thisGap = gap + sigma*d(gen);
        double thisOmega = abs(omega + sigma*d(gen));

        if (thisPlace<M_PI)
            thisPlace = M_PI;
        if (thisGap<0.1)
            thisGap = 0.1;
        if (thisGap>M_PI/2)
            thisGap = M_PI/2;


        double val = funWithGap(thisOmega, thisPlace, thisGap, N);

        if (val < bestVal) {
            success++;
            bestVal = val;
            place = thisPlace;
            gap = thisGap;
            omega = thisOmega;
        }

        if (i%20 == 0) {
            if (success < 5)
                sigma *= c;
            else
                sigma /= c;
            success = 0;
        }
        if (i%5000==0 && (i-1)%5000!=0) {
            cout << "* " << i << ", " << i << ", cv2: " << bestVal << ", sigma: " << sigma << endl;
            cout << "* place = " << place << endl;
            cout << "* gap = " << gap << endl;
            cout << "* omega = " << omega << endl;
        }
    }
    cout << "Initial cv2(" << N << ") = " << bestVal << endl;

    double delta = (2.0*M_PI-gap) / N;
    int i = floor(place/delta+0.5);
    vector<mpreal> phi;
    for (int j=1; j<=i; j++)
        phi.push_back((j-0.5)*delta-M_PI);
    for (int j=i+1; j<=N; j++)
        phi.push_back((j-0.5)*delta-M_PI+gap);

    return make_pair(omega, phi);
}

static int globN;

double gap_cma_objective(const double* par) {
    mpreal::set_default_prec(mpfr::digits2bits(precision));
    double omega = par[0];
    double place = fabs(par[1]);
    double gap = fabs(par[2]);
    while (place<M_PI)
        place += M_PI;
    while (gap>M_PI/2)
        gap -= M_PI/2;
    while (place+gap>2*M_PI)
        place -= M_PI/2;

    double delta = (2.0*M_PI-gap) / globN;
    int i = floor(place/delta+0.5);
    vector<mpreal> params;
    params.push_back(omega);
    for (int j=1; j<=i; j++)
        params.push_back((j-0.5)*delta-M_PI);
    for (int j=i+1; j<=globN; j++)
        params.push_back((j-0.5)*delta+gap-M_PI);
    return cv2(params).toDouble();
}

mpreal cma_objective_mp(const mpreal* par) {
    vector<mpreal> params(globN+1);
    for (int i=0; i<globN+1; i++)
        params[i] = par[i];
    return cv2(params);
}

double cma_objective(const double* par) {
    mpreal::set_default_prec(mpfr::digits2bits(precision));
    vector<mpreal> params(globN+1);
    for (int i=0; i<globN+1; i++)
        params[i] = mpreal(par[i]);
    return cv2(params).toDouble();
}

static FitFunc bipop_cma_objective = [](const double *par, const int N) {
    return cma_objective(par);
};

void gap_cmaes() {

    double init[3] = {0.5, 4.0, 2.0};

    CMA<double,double> cma(3, gap_cma_objective, 3.0, init, false, true);
    double old = cma.mean_y;
    for (int i=0; i<20000; i++) {
//        cout << "Step " << i << endl;
        cma.step();
        if (i%20==0) {
            if (abs(cma.mean_y-old)/old<1e-8)
                break;
            old = cma.mean_y;
//            cout << cma << endl;
        }
    }

    // format output
    cout << "cv2O" << globN << "=" << cma.mean_y << endl;
    cout << "omegaO" << globN << "=" << fabs(cma.mean[0]) << endl;
    cout << "phiO" << globN << "={";

    double place = fabs(cma.mean[1]);
    double gap = fabs(cma.mean[2]);
    while (place<M_PI)
        place += M_PI;
    while (gap>M_PI/2)
        gap -= M_PI/2;
    while (place+gap>2*M_PI)
        place -= M_PI/2;
    double delta = (2.0*M_PI-gap) / globN;
    int i = floor(place/delta+0.5);
    vector<double> phi;
    for (int j=1; j<=i; j++)
        phi.push_back((j-0.5)*delta-M_PI);
    for (int j=i+1; j<=globN; j++)
        phi.push_back((j-0.5)*delta+gap-M_PI);

    for (int i=0; i<globN; i++) {
        double val = phi[i];
        if (val<0)
            val -= 2*M_PI*(floor(-val/(2*M_PI))-1);
        else if (val>2*M_PI)
            val -= 2*M_PI*floor(val/(2*M_PI));
        phi[i] = val;
    }
    sort(phi.begin(), phi.end());
    for (int i=0; i<globN; i++) {
        cout << phi[i];
        if (i==globN-1)
            cout << "}" << endl;
        else
            cout << ", ";
    }
 }

void standard_cmaes() {

    double init[globN+1];
    init[0] = 0.5;
    for (int i=1; i<globN+1; i++)
        init[i] = 2.0*M_PI*rand()/RAND_MAX;

    CMA<double,double> cma(globN+1, cma_objective, 3.0, init, false, true);
    double old = cma.mean_y;
    for (int i=0; i<20000; i++) {
        cma.step();        
        if (i%20==0) {
            if (abs(cma.mean_y-old)/old<1e-8)
                break;
            old = cma.mean_y;
        }
    }

    // format output
    cout << "cv2O" << globN << "=" << cma.mean_y << endl;
    cout << "omegaO" << globN << "=" << cma.mean[0] << endl;
    cout << "phiO" << globN << "={";
    vector<double> phi;
    for (int i=0; i<globN; i++) {
        double val = cma.mean[i+1];
        if (val<0)
            val -= 2*M_PI*(floor(-val/(2*M_PI))-1);
        else if (val>2*M_PI)
            val -= 2*M_PI*floor(val/(2*M_PI));
        phi.push_back(val);
    }
    sort(phi.begin(), phi.end());
    for (int i=0; i<globN; i++) {
        cout << phi[i];
        if (i==globN-1)
            cout << "}" << endl;
        else
            cout << ", ";
    }
}

ProgressFunc<CMAParameters<>,CMASolutions> select_time = [](const CMAParameters<> &cmaparams, const CMASolutions &cmasols)
{
  if (cmasols.niter() % 100 == 0)
    std::cerr << cmasols.elapsed_last_iter() << ": " << cmasols << std::endl;
  return 0;
};

void bipop_cmaes() {

    std::vector<double> init(globN+1);
    init[0] = 0.5;
    for (int i=1; i<globN+1; i++)
        init[i] = 2.0*M_PI*rand()/RAND_MAX;

    double pop_size_multiplier = 3.0;
    CMAParameters<> cmaparams(init,0.1,pop_size_multiplier*(4+3*log((double)(globN+1))));
    cmaparams.set_restarts(5);
    cmaparams.set_algo(BIPOP_CMAES);
//    CMASolutions cmasols = cmaes<>(bipop_cma_objective,cmaparams,select_time);
    CMASolutions cmasols = cmaes<>(bipop_cma_objective,cmaparams);
    std::cout << "* " << cmasols << std::endl;
    std::cout << "* optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n";

    // format output
    cout << "cv2O" << globN << "=" << cmasols.best_candidate().get_fvalue() << endl;
    cout << "omegaO" << globN << "=" << fabs(cmasols.best_candidate().get_x()[0]) << endl;
    cout << "phiO" << globN << "={";
    vector<double> phi = cmasols.best_candidate().get_x();
    phi.erase(phi.begin());
    for (int i=0; i<phi.size(); i++)
        if (phi[i]<0)
            phi[i] -= 2*M_PI*(floor(-phi[i]/(2*M_PI))-1);
        else if (phi[i]>2*M_PI)
            phi[i] -= 2*M_PI*floor(phi[i]/(2*M_PI));
    sort(phi.begin(), phi.end());
    for (int i=0; i<globN; i++) {
        cout << phi[i];
        if (i==globN-1)
            cout << "}" << endl;
        else
            cout << ", ";
    }
}

int main()
{
    srand((int)time(NULL));

      for (int N=10; N<=400; N+=10) {
        precision = 50+ceil(1.487+N*0.647);
        mpreal::set_default_prec(mpfr::digits2bits(precision));
        globN = N;

//      Select appropriate optimization method: bipop is the slowest but gives the best results, gap_cmaes is the fast heuristic proposed in our paper
//        bipop_cmaes();
//        standard_cmaes();
        gap_cmaes();
    }

    return 0;
}

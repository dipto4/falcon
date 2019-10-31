int hermite(real x[NMAX][NDIM], real xdot[NMAX][NDIM], real f[NMAX][NDIM], real fdot[NMAX][NDIM],
        real step[ ], real tlast[ ], real m[ ], real *t, real tend, int n)
{
    real xtemp[NMAX][NDIM], xdottemp[NMAX][NDIM], fi[NDIM],
            fidot[NDIM];
    real a2[NDIM], a3[NDIM], modfi, modfidot, moda2, moda3; real tmin, dt, dts, dtc, dnmtr, temp;
    int imin, nsteps;
    assert(n<=NMAX);
    nsteps=0;
    do {
        tmin = pow(10,6);
        imin = 0;
        loop(i, n) {
            if (step[i] + tlast[i] < tmin) { tmin = step[i] + tlast[i];
            imin = i;
            } 
        }
        assert(imin>=0); *t = tmin;
        i = imin; loop(j, n) {
            dt = *t - tlast[j];
            dts = pow(dt, 2);
            dtc = dts*dt;
            loop(k, NDIM) {
                xtemp[j][k] = x[j][k] + xdot[j][k]*dt + f[j][k]*dts/2 + fdot[j][k]*dtc/6;
                xdottemp[j][k] = xdot[j][k] + f[j][k]*dt + fdot[j][k]*dts/2;
                
            } 
        }
        ffdot(xtemp,xdottemp,m,i,n,fi,fidot); modfi = modfidot = moda2 = moda3 = 0; assert(step[i]!=0);
        loop(k, NDIM) {
            a2[k] = (-6*(f[i][k] - fi[k]) - step[i]*(4*fdot[i][k] + 2*fidot[k]))/pow(step[i],2);
            a3[k] = (12*(f[i][k] - fi[k]) + 6*step[i]*(fdot[i][k] + fidot[k]))/pow(step[i],3);
            x[i][k] = xtemp[i][k] + pow(step[i],4)*a2[k]/24
                + pow(step[i],5)*a3[k]/120;
            xdot[i][k] = xdottemp[i][k] + pow(step[i],3)*a2[k]/6
                + pow(step[i],4)*a3[k]/24;
            f[i][k] = fi[k];
            fdot[i][k] = fidot[k];
            modfi += pow(fi[k],2);
            modfidot += pow(fidot[k],2);
            moda2 += pow(a2[k],2);
            moda3 += pow(a3[k],2);
            
        }
        dnmtr = moda2 + sqrt(modfidot*moda3);
        assert(dnmtr!=0);
        temp = sqrt(0.02*(sqrt(modfi*moda2) + modfidot)/dnmtr); step[i] = min(1.2*step[i], temp);
        tlast[i] = *t;
        nsteps++;
        
    }
    while(*t<tend); return nsteps;
    
}

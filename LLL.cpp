#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;
using ll = long long;
using ld = long double;

class VecR;
class VecZ;

const ld delta = 0.75;

long long round_nearest_away(long double x) 
{
    long long floor_x = static_cast<long long>(floor(x));
    long double diff = x - floor_x;

    if (diff < 0.5L) 
    {
        return floor_x;
    } 
    else if (diff > 0.5L)
    {
        return floor_x + 1;
    } 
    else 
    {
        return (x >= 0.0L) ? floor_x + 1 : floor_x;
    }
}

class VecR
{
private:
    ll m;
    vector<ld> x;
public:
    VecR(){}
    VecR(ll _m): m(_m){x.resize(m);}
    VecR operator-(const VecR& v) const
    {
        VecR v2(m);
        for (ll i = 0; i < m; i++)
            v2.x[i] = x[i] - v.x[i];
        return v2;
    }

    VecR operator+(const VecR& v) const
    {
        VecR v2(m);
        for (ll i = 0; i < m; i++)
            v2.x[i] = x[i] + v.x[i];
        return v2;
    }

    VecR operator/(ld num)
    {
        VecR v(m);
        for (ll i = 0; i < m; i++)
            v.x[i] = x[i]/num;
        return v;
    }

    VecR operator=(const VecZ& v);

    friend VecR operator*(ld d, const VecR& v);
    friend ld inner_product(const VecR& v1, const VecR& v2);
    friend ostream& operator<<(ostream& os, const VecR& v);

    friend ld inner_product(const VecZ& v1, const VecR& v2);
};

VecR operator*(ld d, const VecR& v)
{
    VecR v2(v.m);
    for (ll i = 0; i < v.m; i++)
        v2.x[i] = d * v.x[i];
    return v2;
}

ld inner_product(const VecR& v1, const VecR& v2)
{
    ld d = 0;
    for (ll i = 0; i < v1.m; i++)
        d = d + (v1.x[i] * v2.x[i]);
    return d;
}

ostream& operator<<(ostream& os, const VecR& v)
{
    for (ll i = 0; i < v.m; i++)
        os << v.x[i] << ' ';
    return os;
}


class VecZ
{
private:
    ll m;
    vector<ll> x;
public:
    VecZ(){}
    VecZ(ll _m): m(_m) {x.resize(m);}
    VecZ operator-(const VecZ& v1) const
    {
        VecZ v2(m);
        for (ll i = 0; i < m; i++)
            v2.x[i] = x[i] - v1.x[i];
        return v2;
    }

    friend ostream& operator<<(ostream& os, const VecZ& v);
    friend istream& operator>>(istream& is, VecZ& v);
    friend VecZ operator*(ll c, const VecZ& v);

    friend ld inner_product(const VecZ& v1, const VecR& v2);
    friend __int128 inner_product(const VecZ& v1, const VecZ& v2);

    friend VecR VecR::operator=(const VecZ& v);
};

ostream& operator<<(ostream& os, const VecZ& v)
{
    for (ll i = 0; i < v.m; i++)
        cout << v.x[i] << ' ';
    return os;
}

istream& operator>>(istream& is, VecZ& v)
{
    for (ll i = 0; i < v.m; i++)
        cin >> v.x[i];
    return is;
}

VecZ operator*(ll c, const VecZ& v)
{
    VecZ v2(v.m);
    for (ll i = 0; i < v.m; i++)
        v2.x[i] = c * v.x[i];
    return v2;
}

VecR VecR::operator=(const VecZ& v)
{
    m = v.m;
    for (ll i = 0; i < m; i++)
        x[i] = static_cast<ld>(v.x[i]);
    return *this;
}

__int128 inner_product(const VecZ& v1, const VecZ& v2)
{
    __int128 r = 0;
    for (ll i = 0; i < v1.m; i++)
    {
        __int128 r1 = v1.x[i];
        __int128 r2 = v2.x[i];
        r = r + r1 * r2;
    }
    return r;
}

ld inner_product(const VecZ& v1, const VecR& v2)
{
    ld d = 0;
    for (ll i = 0; i < v1.m; i++)
        d = d + ld(v1.x[i]) * v2.x[i];
    return d;
}

bool isLLLreduced(vector<VecZ>& b, ll n, ll m)
{
   vector<vector<__int128>> gm(n, vector<__int128>(n));
    for (ll i = 0; i < n; i++)
    {
        for (ll j = 0; j <= i; j++)
            gm[i][j] = inner_product(b[i], b[j]);
    }
    vector<vector<ld>> mu(n, vector<ld>(n));
    vector<vector<ld>> r(n, vector<ld>(n));
    vector<vector<ld>> s(n, vector<ld>(n));
    for (ll i = 0; i < n; i++)
    {
        for (ll j = 0; j < i; j++)
        {
             r[i][j] = gm[i][j];
             for (ll k = 0; k < j; k++)
                r[i][j] = r[i][j] - mu[j][k] * r[i][k];
            mu[i][j] = r[i][j]/r[j][j];
        }
        s[i][0] = gm[i][i];
        for (ll j = 1; j <= i; j++)
            s[i][j] = s[i][j - 1] - mu[i][j - 1] * r[i][j - 1];
        r[i][i] = s[i][i];
    }

    for (ll i = 0; i < n; i++)
    {
        for (ll j = 0; j < i; j++)
        {
            if (abs(mu[i][j]) > 0.5L)
                return false;
        }
    }

    for (ll i = 1; i < n; i++)
    {
        if (delta * r[i - 1][i - 1] > s[i][i - 1])
            return false;
    }
    return true;
}

int main()
{

    ll m, n; // Dimension = m, Rank = n
    //cout << "Enter the Dimension and Rank of the Lattice: ";
    cin >> m >> n;

    vector<VecZ> b(n, m); // Basis
    //cout << "Enter the basis:\n";
    for (ll i = 0; i < n; i++)
        cin >> b[i];
    
    /*=====Calculating Gram-Schmidt Vectors=====*/ 
    vector<VecR> _b(n, m); // Gram-Schmidt Vectors
    vector<vector<ld>> mu(n, vector<ld>(n, 0)); // Gram-Schmidt Coefficients
    vector<ld> B(n); // B[i] = Square of Norm of Gram-Schmidt Vector _b[i]
    for (ll i = 0; i < n; i++)
        mu[i][i] = 1;
    _b[0] = b[0];
    B[0] = inner_product(_b[0], _b[0]);
    for (ll i = 1; i < n; i++)
    {
        _b[i] = b[i];
        for (ll j = 0; j < i; j++)
        {
            mu[i][j] = inner_product(b[i], _b[j])/B[j];
            _b[i] = _b[i] - mu[i][j] * _b[j];
        }
        B[i] = inner_product(_b[i], _b[i]);
    }
    /*========LLL Algorithm========*/
    ll ct = 0; // iteration count
    ll k = 1;
    while (k < n)
    {
        //Perform size reduction
        for (ll j = k - 1; j >= 0; j--)
        {
            if (abs(mu[k][j]) > 0.5)
            {
                ll q = round_nearest_away(mu[k][j]);
                b[k] = b[k] - q * b[j];
                mu[k][j] = mu[k][j] - q;
                for (ll i = 0; i < j; i++)
                    mu[k][i] = mu[k][i] - q * mu[j][i];
            }
        }
        // Check for Lovasz condition
        if (B[k] >= (delta - mu[k][k - 1] * mu[k][k - 1]) * B[k - 1])
            k++;
        else
        {
            swap(b[k], b[k - 1]);

            VecR old1 = _b[k - 1];
            VecR old2 = _b[k];
            ld oldB1 = B[k - 1];
            ld oldB2 = B[k];
            _b[k - 1] = old2 + mu[k][k - 1] * old1;
            B[k - 1] = oldB2 + mu[k][k - 1] * mu[k][k - 1] * oldB1;
            _b[k] = ((oldB2 * old1)- ((mu[k][k - 1] * oldB1) * old2))/B[k - 1];
            B[k] = (oldB1 * oldB2)/B[k - 1];

            for (ll j = 0; j < k; j++)
                mu[k][j] = inner_product(b[k], _b[j])/B[j];
            for (ll j = 0; j < k - 1; j++)
                mu[k - 1][j] = inner_product(b[k - 1], _b[j])/B[j];
            
            for (ll i = k + 1; i < n; i++)
                mu[i][k] = inner_product(b[i], _b[k])/B[k];
            for (ll i = k; i < n; i++)
                mu[i][k - 1] = inner_product(b[i], _b[k - 1])/B[k - 1];
            
            k--;
            if (k < 1)
                k = 1;
        } 
        ct++;
        if (!isLLLreduced(b, k, m))
        {
            cerr << "Failure" << endl;
            cerr << ct << endl;
            exit(1);
        }
    }

    if (!isLLLreduced(b, n, m))
    {
        cerr << "Failure" << endl;
        //cerr << ct << endl;
        exit(1);
    }

    //cout << "LLL - reduced basis is:" << '\n';
    for (ll i = 0; i < n; i++)
        cout << b[i] << '\n';
}   

#include "fluidsolver.h"
#include <cmath>

#include <iostream>

using namespace std;

FluidSolver::FluidSolver(int n, float dt)
{
    visc = 0.0f;
    diff = 0.0f;
    this->n = n;
    this->dt = dt;
    size = (n + 2) * (n + 2);
    initArrays();
}

FluidSolver::~FluidSolver()
{
    delete [] u;
    delete [] v;
    delete [] uOld;
    delete [] vOld;
    delete [] d;
    delete [] dOld;
}

float FluidSolver::computeCurl(int i, int j)
{
    float du_dy = (u[I(i, j + 1)] - u[I(i, j - 1)]) * 0.5f;
    float dv_dx = (v[I(i + 1, j)] - v[I(i - 1, j)]) * 0.5f;

    return du_dy - dv_dx;
}

void FluidSolver::vorticityConfinement(float* Fvc_x, float* Fvc_y)
{
    float dw_dx, dw_dy;
    float length;
    float v;

    float vorticityConstant = 1.5f;

    // Calculate magnitude of curl(u,v) for each cell. (|w|)
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            curl[I(i, j)] = std::abs(computeCurl(i, j));
        }
    }

    for (int i = 2; i < n; i++)
    {
        for (int j = 2; j < n; j++)
        {
            // Find derivative of the magnitude (n = del |w|)
            dw_dx = (curl[I(i + 1, j)] - curl[I(i - 1, j)]) * 0.5f;
            dw_dy = (curl[I(i, j + 1)] - curl[I(i, j - 1)]) * 0.5f;

            // Calculate vector length. (|n|)
            // Add small factor to prevent divide by zeros.
            length = (float) std::sqrt(dw_dx * dw_dx + dw_dy * dw_dy)
                     + 0.000001f;

            // N = ( n/|n| )
            dw_dx /= length;
            dw_dy /= length;

            v = computeCurl(i, j);

            // N x w
            Fvc_x[I(i, j)] = vorticityConstant * dw_dy * -v;
            Fvc_y[I(i, j)] = vorticityConstant * dw_dx *  v;
        }
    }
}

void FluidSolver::setBoundary(int b, float* x)
{
    for (int i = 1; i <= n; i++)
    {
        x[I(  0, i  )] = b == 1 ? -x[I(1, i)] : x[I(1, i)];
        x[I(n+1, i  )] = b == 1 ? -x[I(n, i)] : x[I(n, i)];
        x[I(  i, 0  )] = b == 2 ? -x[I(i, 1)] : x[I(i, 1)];
        x[I(  i, n+1)] = b == 2 ? -x[I(i, n)] : x[I(i, n)];
    }

    x[I(  0,   0)] = 0.5f * (x[I(1, 0  )] + x[I(  0, 1)]);
    x[I(  0, n+1)] = 0.5f * (x[I(1, n+1)] + x[I(  0, n)]);
    x[I(n+1,   0)] = 0.5f * (x[I(n, 0  )] + x[I(n+1, 1)]);
    x[I(n+1, n+1)] = 0.5f * (x[I(n, n+1)] + x[I(n+1, n)]);

}

void FluidSolver::advect(int b, float* d, float* d0, float* du, float* dv)
{
    int i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt * n;

    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            // go backwards through velocity field
            x = i - dt0 * du[I(i, j)];
            y = j - dt0 * dv[I(i, j)];

            // interpolate results
            // first, clamp x
            if (x > n + 0.5) x = n + 0.5f;
            if (x < 0.5)     x = 0.5f;

            // find the cell x is in and the next cell
            i0 = (int) x;
            i1 = i0 + 1;

            // similarly, clamp y
            if (y > n + 0.5) y = n + 0.5f;
            if (y < 0.5)     y = 0.5f;

            // find the cell y is in and the next cell
            j0 = (int) y;
            j1 = j0 + 1;

            // now do the interpolation
            // first, find the interpolation coefficients
            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;

            // bilinear interpolation to update the new velocity at cell(i,j)
            // in each term, take the two components away from the point whose value you're using
            d[I(i, j)] = s0 * (t0 * d0[I(i0, j0)] + t1 * d0[I(i0, j1)])
                       + s1 * (t0 * d0[I(i1, j0)] + t1 * d0[I(i1, j1)]);

        }
    }
    setBoundary(b, d);
}

void FluidSolver::diffuse(int b, float* c, float* c0, float diff)
{
    // called in velocitySolver as diffuse(0, u, uOld, visc); diffuse(0, v, vOld, visc);
    float a = dt * diff * n * n;
    linearSolver(b, c, c0, a, 1 + 4 * a);
}

void FluidSolver::swap(float* a, float* b)
{
    for (int i = 0; i < size; i++)
    {
        float tmp = a[i];
        a[i] = b[i];
        b[i] = tmp;
    }
}

void FluidSolver::addSource(float* a, float* b)
{
    for (int i = 0; i < size; i++)
        a[i] += dt * b[i];
}

void FluidSolver::initArrays()
{
    d    = new float[size];
    dOld = new float[size];
    u    = new float[size];
    uOld = new float[size];
    v    = new float[size];
    vOld = new float[size];
    curl = new float[size];

    for (int i = 0; i < size; i++)
    {
        u[i] = uOld[i] = v[i] = vOld[i] = 0.0f;
        d[i] = dOld[i] = curl[i] = 0.0f;
    }
}

void FluidSolver::buoyancy(float* fbuoy)
{
    float Tamb = 0;
    float a = 0.00250f;
    // float a = 0.000625f;
    float b = 0.025f;

    // sum all temperatures
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            Tamb += d[I(i, j)];
        }
    }

    // get average temperature
    Tamb /= (n * n);

    // for each cell compute buoyancy force
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            // fbuoy[I(i, j)] = a * d[I(i, j)] + -b * (d[I(i, j)] - Tamb);
            fbuoy[I(i, j)] = a * d[I(i, j)];
        }
    }
}

void FluidSolver::step()
{
    temperatureSolver();
    velocitySolver();
    densitySolver();
    implicitSurfaceSolver();
}

void FluidSolver::temperatureSolver()
{}

void FluidSolver::implicitSurfaceSolver()
{}

// page 5, stam - stable fluids
void FluidSolver::velocitySolver()
{
    // cout << "velocity solver " << endl;
    vorticityConfinement(uOld, vOld);
    addSource(u, uOld);
    addSource(v, vOld);

    buoyancy(vOld);
    addSource(v, vOld);

    swap(u, uOld);
    diffuse(0, u, uOld, visc);
    swap(v, vOld);
    diffuse(0, v, vOld, visc);

    project(u, v, uOld, vOld);
    
    swap(u, uOld);
    swap(v, vOld);

    advect(1, u, uOld, uOld, vOld);
    advect(2, v, vOld, uOld, vOld);

    project(u, v, uOld, vOld);

    clear(uOld); clear(vOld);
}

void FluidSolver::project(float* x, float* y, float* p, float* div)
{
    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            // division by 2 is done because we are doing central differencing
            div[I(i, j)] = (x[I(i+1, j)] - x[I(i-1, j)]
                          + y[I(i, j+1)] - y[I(i, j-1)])
                          * - 0.5f / n;
            p[I(i, j)] = 0;
        }
    }

    setBoundary(0, div);
    setBoundary(0, p);

    // now obtain the pressure by running the linear solver
    linearSolver(0, p, div, 1, 4);

    for (int i = 1; i <= n; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            x[I(i, j)] -= 0.5f * n * (p[I(i+1, j)] - p[I(i-1, j)]);
            y[I(i, j)] -= 0.5f * n * (p[I(i, j+1)] - p[I(i, j-1)]);
        }
    }

    setBoundary(1, x);
    setBoundary(2, y);
}

void FluidSolver::linearSolver(int b, float* x, float* x0, float a, float c)
{
    for (int k = 0; k < 20; k++)
    {
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n; j++)
            {
                x[I(i, j)] = (a * ( x[I(i-1, j)] + x[I(i+1, j)]
                                +   x[I(i, j-1)] + x[I(i, j+1)])
                                +  x0[I(i, j)]) / c;
            }
        }
        setBoundary(b, x);
    }
}

void FluidSolver::initFluidSingle(pair<int,int> location)
{
    int x = location.first;
    int y = location.second;
    d[I(x,y)] += densityIncrement;
}

void FluidSolver::initFluid(pair<int,int> locations[], int num_locations)
{
    for (int i = 0 ; i < num_locations ; i++)
    {
        initFluidSingle(locations[i]);
    }
}

void FluidSolver::densitySolver()
{
    // cout << "density solver " << endl;
    addSource(d, dOld);
    swap(d, dOld);
    diffuse(0, d, dOld, diff);
    swap(d, dOld);
    advect(0, d, dOld, u, v);
    clear(dOld);
}
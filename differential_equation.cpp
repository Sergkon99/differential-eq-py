#include <stdexcept>
#include <iostream>
#include <set>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <sstream>

#include "differential_equation.hpp"

DifferentialEquation::DifferentialEquation(FuncType init_func, double h, int T) :
    _init_func(init_func), _h(h), _T(T)
{
    for (double t = -_T; t <= 0; t += h)
    {
        _x.push_back(_init_func(t));
        _t.push_back(t);
    }
    _offset = _x.size();
}

void DifferentialEquation::Add(double val, double t)
{
    _x.push_back(val);
    _t.push_back(t);
}

double DifferentialEquation::Call(double val)
{
    return _init_func(val);
}

double DifferentialEquation::GetCur()
{
    return _x.back();
}

double DifferentialEquation::GetWithDelay()
{
    if (_x.size() >= _offset)
        return _x[_x.size() - _offset];
    else
        throw std::runtime_error("Index out of range for delay");
}

Coords DifferentialEquation::GetSolutionData()
{
    return { _t, _x };
}

size_t DifferentialEquation::GetSize() const
{
    return _x.size();
}

double DifferentialEquation::Get(size_t i) const
{
    return _x[i];
}


DifferentialEquationSystem::DifferentialEquationSystem(FuncList init_funcs, double h, int T, int L, RelationMatrix relations, FuncType finit_func) :
    _h(h), _T(T), _L(L), _t(0), _finit_func(finit_func)
{
    _relation_map = RelationBuilder(relations).Get();
    _equations.clear();
    for (FuncType init_func : init_funcs)
        _equations.push_back(DifferentialEquation(init_func, _h, _T));
    _n = _equations.size();
    if (relations.empty() || _n != relations.size() || _n != relations[0].size())
    {
        std::stringstream ss;
        ss << "Incorrect relations size. Need " << _n << ". Found " << relations.size();
        throw std::runtime_error(ss.str());
    }
}

CoordsList DifferentialEquationSystem::GetSolutionData()
{
    CoordsList result;
    for (auto e : _equations)
        result.push_back(e.GetSolutionData());
    return result;
}

double DifferentialEquationSystem::F(double cur, double delay, size_t i, bool use_force)
{
    double ans = -cur + _L * _finit_func(delay);
    if (use_force)
        ans += D(i);
    return ans;
}

double DifferentialEquationSystem::D(size_t i) 
{
    // i - порядковый номер решения
    double _bond_force = 0;
    if (_relation_map.find(i) == _relation_map.end())
        return 0;

    for (auto edge : _relation_map.at(i)) 
    {
        int num = edge.first;
        double w = edge.second;
        _bond_force += w * (_equations[num].GetCur() - _equations[i].GetCur());
    }
    return _bond_force;
}

void DifferentialEquationSystem::Solve(double end_time) 
{
    int solve_per_cent = 0;
    while (_t < end_time) 
    {
        std::vector< double > next_xs(_n);
        for (size_t i = 0; i < _n; i++) 
        {
            double cur_x = _equations[i].GetCur(),
                delay_x = _equations[i].GetWithDelay(),
                k1 = F(cur_x, delay_x, i),
                k2 = F(cur_x + _h / 2 * k1, delay_x, i),
                k3 = F(cur_x + 1 / 2 * k2, delay_x, i),
                k4 = F(cur_x + 1 * k3, delay_x, i)
                ;
            double next_x = cur_x + _h * k2;
            next_xs[i] = next_x;
        }

        _t += _h;
        std::vector< double > tmp;
        for (size_t i = 0; i < _n; i++) 
        {
            _equations[i].Add(next_xs[i], _t);
        }
    }
}
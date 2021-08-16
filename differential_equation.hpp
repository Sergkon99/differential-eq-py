#pragma once
#include <vector>
#include <map>
#include <deque>
#include <cassert>
#include <optional>
#include <functional>

typedef std::pair< std::vector< double >, std::vector< double > > Coords;
typedef std::vector< Coords > CoordsList;

typedef std::map< size_t, std::vector< std::pair< int, double > > > RelationMap;
typedef std::vector< std::vector< double > > RelationMatrix;
typedef std::function<double(double)> FuncType;
typedef std::vector< FuncType > FuncList;


class RelationBuilder
{
    // ������� ������
    RelationMatrix _matrix;
    // ������� ������
    RelationMap _data;
    // ����� ���������
    size_t _count;
public:
    RelationBuilder(RelationMatrix matrix) :
        _matrix(matrix)
    {
        assert(!_matrix.empty() && _matrix[0].size() == _matrix.size() && "Incorrect relation size");
        _count = matrix.size();
        for (size_t i = 0; i < _count; i++)
        {
            for (size_t j = 0; j < _count; j++)
            {
                if (_matrix[i][j] != 0.0)
                {
                    _data[i].push_back({ j, _matrix[i][j] });
                }
            }
        }
    };
    RelationMap Get()
    {
        return _data;
    }
};

class DifferentialEquation
{
private:
    // ��������� �������(������� �������������)
    FuncType _init_func;
    // �������� ��������� - ��� ��� �����
    double _h;
    // �������� ��������� - ������������
    int _T;
    // ������ ����� ���������� �������
    std::vector< double > _x, _t;
    // ���-�� ����� � ������������� ��� ��������� �����
    size_t _offset;
public:
    DifferentialEquation(FuncType init_func, double h, int T);

    /*
    * \brief �������� ����� � ���������
    * \param val ��������
    * \param t �����
    **/
    void Add(double val, double t);
    /*
    * \brief ������� ������� �������������
    * \param val ��������
    **/
    double Call(double val);
    /*
    * \brief ������� ������� �������� ��������
    **/
    double GetCur();
    /*
    * \brief ������� �������� �������� � ������ ������������
    **/
    double GetWithDelay();
    /*
    * \brief ������� �������� ����������� �������
    */
    Coords GetSolutionData();
    /*
    * \brief ������� ����� ����� � �������
    */
    size_t GetSize() const;
    /*
    * \brief ������� �������� ������� �� �������
    */
    double Get(size_t i) const;
};

class DifferentialEquationSystem
{
private:
    // 
    double EPS = 0.00001;
    // ��������� �������
    std::vector< DifferentialEquation > _equations;
    // �������� ��������� - ��� ��� �����
    double _h;
    // �������� ��������� - ������������
    int _T;
    // �������� ���������
    int _L;
    // �����
    double _t;
    // ���-�� ������������ � �������
    size_t _n;
    // ������ ������
    RelationMap _relation_map;
    // �������� ��������� - ������� F
    FuncType _finit_func;

    // ������� F
    double F(double cur, double delay, size_t i, bool use_force = true);
    // ���� �����
    double D(size_t i);
    double D(double cur, double delay, size_t i);
public:
    DifferentialEquationSystem(std::vector< FuncType >, double h, int T, int L, RelationMatrix relations, FuncType finit_func);

    /*
    * \brief ������ �������
    * \param end_time �� ������ ������� ������ �������
    **/
    void Solve(double end_time);
    /*
    * \brief ������� �������� ���� ���������� �������
    */
    CoordsList GetSolutionData();
};

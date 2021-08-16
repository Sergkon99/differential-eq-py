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
    // Матрица связей
    RelationMatrix _matrix;
    // Словарь связей
    RelationMap _data;
    // Число уравнений
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
    // Начальная функция(функция инициализации)
    FuncType _init_func;
    // Параметр уравнения - шаг при счете
    double _h;
    // Параметр уравнения - запаздывание
    int _T;
    // Список точек численного решения
    std::vector< double > _x, _t;
    // Кол-во точек в запаздываении при численном счете
    size_t _offset;
public:
    DifferentialEquation(FuncType init_func, double h, int T);

    /*
    * \brief Добавить точку в уравнение
    * \param val Значение
    * \param t Время
    **/
    void Add(double val, double t);
    /*
    * \brief Вызвать функцию инициализации
    * \param val Параметр
    **/
    double Call(double val);
    /*
    * \brief Вернуть текущее значение уравения
    **/
    double GetCur();
    /*
    * \brief Вернуть значение уравения с учетом запаздывания
    **/
    double GetWithDelay();
    /*
    * \brief Вернуть значения сосчитаного решения
    */
    Coords GetSolutionData();
    /*
    * \brief Вернуть число точек в решении
    */
    size_t GetSize() const;
    /*
    * \brief Вернуть значение решения по индексу
    */
    double Get(size_t i) const;
};

class DifferentialEquationSystem
{
private:
    // 
    double EPS = 0.00001;
    // Численные решения
    std::vector< DifferentialEquation > _equations;
    // Параметр уравнения - шаг при счете
    double _h;
    // Параметр уравнения - запаздывание
    int _T;
    // Параметр уравнения
    int _L;
    // Время
    double _t;
    // Кол-во осцилляторов в системе
    size_t _n;
    // Список связей
    RelationMap _relation_map;
    // Параметр уравнения - функция F
    FuncType _finit_func;

    // Функция F
    double F(double cur, double delay, size_t i, bool use_force = true);
    // Сила свзяи
    double D(size_t i);
    double D(double cur, double delay, size_t i);
public:
    DifferentialEquationSystem(std::vector< FuncType >, double h, int T, int L, RelationMatrix relations, FuncType finit_func);

    /*
    * \brief Решить систему
    * \param end_time До какого времени решать систему
    **/
    void Solve(double end_time);
    /*
    * \brief Вернуть значения всех сосчитаных решений
    */
    CoordsList GetSolutionData();
};

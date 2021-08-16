#include <iostream>
#include <string>
#include <vector>

#include <boost/python.hpp>

#include "differential_equation.hpp"

using namespace boost::python;

std::function<double(double)> CreateCallback(object callback) 
{
	if (!PyCallable_Check(callback.ptr()))
		throw std::runtime_error("Callback object is not callable");

	auto lambda_func = [callback](double t) 
	{
		return extract<double>(callback(t)); 
	};
	return static_cast<std::function<double(double)>>(lambda_func);
}

std::shared_ptr<DifferentialEquation> DifferentialEquation_Constructor(object init_func_py, double h, int T) 
{
	return std::shared_ptr<DifferentialEquation>(new DifferentialEquation(CreateCallback(init_func_py), h, T ));
}

list DifferentialEquation_GetSolutionData(DifferentialEquation& self)
{
	Coords data = self.GetSolutionData();
	list py_list;
	auto xs = data.second, ts = data.first;
	for (size_t i = 0; i < xs.size(); i++)
		py_list.append(make_tuple(ts[i], xs[i]));

	return py_list;
}

std::shared_ptr<DifferentialEquationSystem> DifferentialEquationSystem_Constructor(list init_funcs_py, double h, int T, int L, list relations_py, object finit_func_py) 
{
	std::vector<FuncType> init_funcs;
	for (size_t i = 0; i < len(init_funcs_py); i++) 
	{
		object py_func = extract<object>(init_funcs_py[i]);
		init_funcs.push_back(CreateCallback(py_func));
	}
	RelationMatrix relations;
	for (size_t i = 0; i < len(relations_py); i++) 
	{
		list inside_list = extract<list>(relations_py[i]);
		std::vector<double> rel_inside;
		for (size_t j = 0; j < len(inside_list); j++) 
		{
			rel_inside.push_back(extract<double>(inside_list[j]));
		}
		relations.push_back(rel_inside);
	}
	return std::shared_ptr<DifferentialEquationSystem>(new DifferentialEquationSystem(init_funcs, h, T, L, relations, CreateCallback(finit_func_py)));
}

list DifferentialEquationSystem_GetSolutionData(DifferentialEquationSystem& self)
{
	CoordsList data = self.GetSolutionData();
	list py_list;
	for (auto item : data) 
	{
		list py_list_inside_t, py_list_inside_x;
		auto xs = item.second, ts = item.first;
		for (size_t i = 0; i < xs.size(); i++) 
		{
			py_list_inside_t.append(ts[i]);
			py_list_inside_x.append(xs[i]);
		}
		py_list.append(make_tuple(py_list_inside_t, py_list_inside_x));
	}
	return py_list;
}

BOOST_PYTHON_MODULE(differential_eq_py)
{
	class_<std::shared_ptr<DifferentialEquation>>("DifferentialEquation", no_init)
		.def("__init__", make_constructor(&DifferentialEquation_Constructor, default_call_policies(), (arg_("init_func"), arg_("h"), arg_("t"))))
		.def("GetSolutionData", &DifferentialEquation_GetSolutionData)
		.def("Call", &DifferentialEquation::Call)
		;
	class_<std::shared_ptr<DifferentialEquationSystem>>("DifferentialEquationSystem", no_init)
		.def("__init__", make_constructor(&DifferentialEquationSystem_Constructor, default_call_policies(), (arg_("init_funcs"), arg_("h"), arg_("T"), arg_("L"), arg_("relations"), arg_("finit_func"))))
		.def("GetSolutionData", &DifferentialEquationSystem_GetSolutionData)
		.def("Solve", &DifferentialEquationSystem::Solve)
		;
}
//#include "opaque.h"
//#include <pybind11/pybind11.h>


//namespace py = pybind11;
//using namespace py::literals;
//
//
//PYBIND11_EMBEDDED_MODULE(opaque_module, m)
//{
//py::bind_vector<std::vector<double>>(m, "VectorDouble");
//}

//PYBIND11_MODULE(opaque, m){
//py::class_<opaque_vector>(m, "DoubleVector")
//    .def(py::init<>())
//    .def("push_back", [](opaque_vector & op, const double &x) {op.V.push_back(x);})
//    .def("__idk__", [](opaque_vector & op, const std::vector<double> &v) {op.V = v;})
//    .def("__getitem__", [](const opaque_vector & op, const std::size_t & idx) {return op.V.at(idx);});
//    .def("clear", &std::vector<double>::clear)
//    .def("pop_back", &std::vector<double>::pop_back)
//    .def("push_back", (void(std::vector<double>::*)(const double &)) & std::vector<double>::push_back);
//    .def("size", (int(std::vector<double>::*)(void)) & std::vector<double>::size)
//    .def("__idk__", [](std::vector<double> &v, std::vector<double> &v_in) {v = v_in;})
//    .def("__getitem__", (double &(std::vector<double>::*)(std::size_t)) & std::vector<double>::at)
//    .def("__len__", [](const std::vector<double> &v) { return v.size(); })
//    .def("__iter__", [](std::vector<double> &v) {
//        return py::make_iterator(v.begin(), v.end());
//    }, py::keep_alive<0, 1>()); /* Keep vector alive while iterator is used */
//}
//
//template <typename T>
//py::array_t<T> wrapper(py::array_t<T> input) {
//    auto proxy = input.template unchecked<1>();
//    std::vector<T> *result = new std::vector<T>(compute_something_returns_vector(proxy));
//
//    py::capsule free_when_done(result, [](void *f) {
//        auto foo = reinterpret_cast<std::vector<T> *>(f);
//        delete foo;
//    });
//
//    return py::array_t<T>({result->size()}, // shape
//                          {sizeof(T)},      // stride
//                          result->data(),   // data pointer
//                          free_when_done);
//}
//
//PYBIND11_MODULE(opaque, m){
//    m.def("my_func", &wrapper<int>);
//    m.def("my_func", &wrapper<double>);
//}
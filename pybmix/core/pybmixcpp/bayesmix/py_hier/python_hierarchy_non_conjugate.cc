#include "python_hierarchy_non_conjugate.h"

#include <google/protobuf/stubs/casts.h>
#include <pybind11/embed.h>
#include <pybind11/pybind11.h>

#include <Eigen/Dense>
#include <random>
#include <sstream>
#include <stan/math/prim/prob.hpp>
#include <string>
#include <vector>
#include <cmath>

#include "algorithm_state.pb.h"
#include "hierarchy_prior.pb.h"
#include "ls_state.pb.h"
#include "bayesmix/src/utils/rng.h"
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "auxiliary_functions.h"

//! PYTHON
/*
double PythonHierarchyNonConjugate::like_lpdf(const Eigen::RowVectorXd &datum) const {
    double result = like_lpdf_evaluator(datum, state.generic_state).cast<double>();
    return result;
}
*/
double PythonHierarchyNonConjugate::like_lpdf(const Eigen::RowVectorXd &datum) const {
  return stan::math::double_exponential_lpdf(datum(0), state.generic_state[0],
                                             state.generic_state[1]);
}


//! PYTHON
//void PythonHierarchyNonConjugate::initialize_state() {
//    py::list state_py = initialize_state_evaluator(hypers->generic_hypers);
//    state.generic_state = list_to_vector(state_py);
//}
void PythonHierarchyNonConjugate::initialize_state() {
  state.generic_state.push_back(hypers->generic_hypers[0]);
  state.generic_state.push_back(hypers->generic_hypers[3] / (hypers->generic_hypers[2] + 1));  // mode of Inv-Gamma
}


//! C++
void PythonHierarchyNonConjugate::initialize_hypers() {
    py::list hypers_py = initialize_hypers_evaluator();
    hypers->generic_hypers = list_to_vector(hypers_py);
}

//! PYTHON
void PythonHierarchyNonConjugate::update_hypers(
        const std::vector <bayesmix::AlgorithmState::ClusterState> &states) {
    auto &rng = bayesmix::Rng::Instance().get();
    if (prior->has_values()) return;
}



//! PYTHON
//Python::State PythonHierarchyNonConjugate::draw(const Python::Hyperparams &params) {
//  Python::State out;
//  synchronize_cpp_to_py_state(bayesmix::Rng::Instance().get(), py_gen);
//  py::list draw_py = draw_evaluator(state.generic_state,params.generic_hypers,py_gen);
//  out.generic_state = list_to_vector(draw_py);
//  synchronize_py_to_cpp_state(bayesmix::Rng::Instance().get(), py_gen);
//  return out;
//}
Python::State PythonHierarchyNonConjugate::draw(const Python::Hyperparams &params) {
  auto &rng = bayesmix::Rng::Instance().get();
  Python::State out{};
  out.generic_state.push_back(stan::math::normal_rng(params.generic_hypers[0], sqrt(params.generic_hypers[1]), rng));
  out.generic_state.push_back(stan::math::inv_gamma_rng(params.generic_hypers[2], 1. / params.generic_hypers[3], rng));
  return out;
}


//! PYTHON
//void PythonHierarchyNonConjugate::update_summary_statistics(
//        const Eigen::RowVectorXd &datum, const bool add) {
//    py::list results = update_summary_statistics_evaluator(datum,add,sum_stats, state.generic_state, cluster_data_values);
//    py::list sum_stats_py = results[0];
//    py::array cluster_data_values_py = results[1];
//    sum_stats = list_to_vector(sum_stats_py);
//    cluster_data_values = cluster_data_values_py.cast<Eigen::MatrixXd>();
//}

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void PythonHierarchyNonConjugate::update_summary_statistics(
    const Eigen::RowVectorXd &datum, const bool add) {
  if (add) {
    sum_stats[0] += std::abs(state.generic_state[0] - datum(0, 0));
    // cluster_data_values.push_back(datum);
    cluster_data_values.conservativeResize(cluster_data_values.rows()+1, cluster_data_values.cols());
    cluster_data_values.col(cluster_data_values.cols()-1) = datum;
  } else {
    sum_stats[0] -= std::abs(state.generic_state[0] - datum(0, 0));
    auto it = std::find(cluster_data_values.begin(), cluster_data_values.end(),
                        datum);
    cluster_data_values.erase(it); //not working
  }
}

//! PYTHON
//void PythonHierarchyNonConjugate::clear_summary_statistics() {
//    py::list sum_stats_py = clear_summary_statistics_evaluator(sum_stats);
//    sum_stats = list_to_vector(sum_stats_py);
//}
void PythonHierarchyNonConjugate::clear_summary_statistics() {
  cluster_data_values.clear();
  sum_stats[0] = 0;
  sum_stats[1] = 0;
}


//! PYTHON
//void PythonHierarchyNonConjugate::sample_full_cond(const bool update_params /*= false*/) {
//  if (this->card == 0) {
//    // No posterior update possible
//    this->sample_prior();
//  } else {
//    synchronize_cpp_to_py_state(bayesmix::Rng::Instance().get(), py_gen);
//    py::list result = sample_full_cond_evaluator(state.generic_state, sum_stats, py_gen, cluster_data_values, hypers->generic_hypers);
//    synchronize_py_to_cpp_state(bayesmix::Rng::Instance().get(), py_gen);
//    py::list state_list = result[0];
//    py::list sum_stats_list = result[1];
//    state.generic_state = list_to_vector(state_list);
//    sum_stats = list_to_vector(sum_stats_list);
//    }
//}
void PythonHierarchyNonConjugate::sample_full_cond(const bool update_params /*= false*/) {
  if (this->card == 0) {
    // No posterior update possible
    this->sample_prior();
  } else {
    // Number of iterations to compute the acceptance rate
    ++iter_;

    // Random generator
    auto &rng = bayesmix::Rng::Instance().get();

    // Candidate mean and candidate log_scale
    Eigen::VectorXd curr_unc_params(2);
    curr_unc_params << state(0), std::log(state(1));

    Eigen::VectorXd prop_unc_params = propose_rwmh(curr_unc_params);

    double log_target_prop =
        eval_prior_lpdf_unconstrained(prop_unc_params) +
        eval_like_lpdf_unconstrained(prop_unc_params, false);

    double log_target_curr =
        eval_prior_lpdf_unconstrained(curr_unc_params) +
        eval_like_lpdf_unconstrained(curr_unc_params, true);

    double log_a_rate = log_target_prop - log_target_curr;

    if (std::log(stan::math::uniform_rng(0, 1, rng)) < log_a_rate) {
      state(0) = prop_unc_params(0);
      state(1) = std::exp(prop_unc_params(1));
      sum_stats(0) = sum_stats(1);
    }
  }
}


//! PYTHON
/*
Eigen::VectorXd PythonHierarchyNonConjugate::propose_rwmh(
    const Eigen::VectorXd &curr_vals) {
    synchronize_cpp_to_py_state(bayesmix::Rng::Instance().get(), py_gen);
    py::list proposal = propose_rwmh_evaluator(curr_vals, hypers->generic_hypers, py_gen);
    synchronize_py_to_cpp_state(bayesmix::Rng::Instance().get(), py_gen);
    double candidate_mean = proposal[0].cast<double>();
    double candidate_log_scale = proposal[1].cast<double>();
    Eigen::VectorXd proposalcpp(2);
    proposalcpp << candidate_mean, candidate_log_scale;
    return proposalcpp;
}
*/
Eigen::VectorXd PythonHierarchyNonConjugate::propose_rwmh(
    const Eigen::VectorXd &curr_vals) {
  auto &rng = bayesmix::Rng::Instance().get();
  double candidate_mean =
      curr_vals(0) + stan::math::normal_rng(0, sqrt(hypers->generic_hypers[4]), rng);
  double candidate_log_scale =
      curr_vals(1) +
      stan::math::normal_rng(0, sqrt(hypers->generic_hypers[5]), rng);
  Eigen::VectorXd proposal(2);
  proposal << candidate_mean, candidate_log_scale;
  return proposal;
}

//! PYTHON
/*
double PythonHierarchyNonConjugate::eval_prior_lpdf_unconstrained(
    const Eigen::VectorXd &unconstrained_parameters) {

    double result = eval_prior_lpdf_unconstrained_evaluator(unconstrained_parameters, hypers->generic_hypers).cast<double>();
    return result;
}
*/
double PythonHierarchyNonConjugate::eval_prior_lpdf_unconstrained(
    const Eigen::VectorXd &unconstrained_parameters) {
  double mu = unconstrained_parameters(0);
  double log_scale = unconstrained_parameters(1);
  double scale = std::exp(log_scale);
  return stan::math::normal_lpdf(mu, hypers->generic_hypers[0], std::sqrt(hypers->generic_hypers[1])) +
         stan::math::inv_gamma_lpdf(scale, hypers->generic_hypers[2], hypers->generic_hypers[3]) +
         log_scale;
}

//! PYTHON
/*
double PythonHierarchyNonConjugate::eval_like_lpdf_unconstrained(
    const Eigen::VectorXd &unconstrained_parameters, const bool is_current) {
    double result = eval_like_lpdf_unconstrained_evaluator(unconstrained_parameters, is_current, sum_stats, cluster_data_values).cast<double>();
    return result;
}
*/

double PythonHierarchyNonConjugate::eval_like_lpdf_unconstrained(
    const Eigen::VectorXd &unconstrained_parameters, const bool is_current) {
  double mean = unconstrained_parameters(0);
  double log_scale = unconstrained_parameters(1);
  double scale = std::exp(log_scale);
  double diff_sum = 0;  // Sum of absolute values of data - candidate_mean
  if (is_current) {
    diff_sum = sum_stats[0];
  } else {
    for (auto &elem : cluster_data_values) {
      diff_sum += std::abs(elem(0, 0) - mean);
    }
    sum_stats[1] = diff_sum;
  }
  return std::log(0.5 / scale) + (-0.5 / scale * diff_sum);
}

//! C++
void PythonHierarchyNonConjugate::set_state_from_proto(
        const google::protobuf::Message &state_) {
    auto &statecast = downcast_state(state_);
    int size = statecast.general_state().size();
    std::vector<double> aux_v{};
    for (int i = 0; i < size; ++i) {
        aux_v.push_back((statecast.general_state().data())[i]);
    }
    state.generic_state = aux_v;
    set_card(statecast.cardinality());
}

//! C++
std::shared_ptr <bayesmix::AlgorithmState::ClusterState>
PythonHierarchyNonConjugate::get_state_proto() const {
    bayesmix::Vector state_;
    state_.set_size(state.generic_state.size());
    *state_.mutable_data() = {
            state.generic_state.data(),
            state.generic_state.data() + state.generic_state.size()};
    auto out = std::make_shared<bayesmix::AlgorithmState::ClusterState>();
    out->mutable_general_state()->CopyFrom(state_);
    return out;
}

//! C++
void PythonHierarchyNonConjugate::set_hypers_from_proto(
        const google::protobuf::Message &hypers_) {
    auto &hyperscast = downcast_hypers(hypers_).python_state();
    int size = hyperscast.data().size();
    std::vector<double> aux_v{};
    for (int i = 0; i < size; ++i) {
        aux_v.push_back((hyperscast.data())[i]);
    }
    hypers->generic_hypers = aux_v;
}

//! C++
std::shared_ptr <bayesmix::AlgorithmState::HierarchyHypers>
PythonHierarchyNonConjugate::get_hypers_proto() const {
    bayesmix::Vector hypers_;
    hypers_.set_size(hypers->generic_hypers.size());
    *hypers_.mutable_data() = {
            hypers->generic_hypers.data(),
            hypers->generic_hypers.data() + hypers->generic_hypers.size()};
    auto out = std::make_shared<bayesmix::AlgorithmState::HierarchyHypers>();
    out->mutable_python_state()->CopyFrom(hypers_);
    return out;
}
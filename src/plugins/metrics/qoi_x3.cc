#include <cmath>
#include "pressio_data.h"
#include "pressio_options.h"
#include "pressio_compressor.h"
#include "libpressio_ext/cpp/data.h"
#include "libpressio_ext/cpp/metrics.h"
#include "libpressio_ext/cpp/options.h"
#include "libpressio_ext/cpp/pressio.h"
#include "libpressio_ext/cpp/domain_manager.h"
#include "std_compat/memory.h"

namespace libpressio {
namespace qoi_x3 {
  struct qoi_x3_metrics {
    double min_val_x3;
    double max_val_x3;
    double max_abs_val_diff_x3;
    double max_rel_val_diff_x3; 
  };

  inline double qoi_fn(double x) {
    return x*x*x;
  }

  struct compute_metrics{
    template <class ForwardIt1, class ForwardIt2>
    qoi_x3_metrics operator()(ForwardIt1 input_begin, ForwardIt1 input_end, ForwardIt2 input2_begin, ForwardIt2 input2_end)
    {
      qoi_x3_metrics m{};
      double max_x3_error = 0; 
      if(input_begin != nullptr && input2_begin != nullptr) {
        double val_x3_min = double(*input_begin) * double(*input_begin); 
        double val_x3_max = val_x3_min; 
        double val_x3_diff_min = std::abs(qoi_fn(double(*input_begin))
                                           -qoi_fn(double(*input2_begin))); 
        double val_x3_diff_max = val_x3_diff_min; 


        while(input_begin != input_end && input2_begin != input2_end) {
          double cur_sq_err = std::abs(qoi_fn(double(*input_begin))
                                         -qoi_fn(double(*input2_begin)));
          double cur_sq = qoi_fn(double(*input_begin));
          // x3 
          val_x3_min = std::min(val_x3_min, cur_sq);
          val_x3_max = std::max(val_x3_max, cur_sq);
          val_x3_diff_min = std::min(val_x3_diff_min, cur_sq_err);
          val_x3_diff_max = std::max(val_x3_diff_max, cur_sq_err);
          ++input_begin;
          ++input2_begin;
        }
        m.max_val_x3 = val_x3_max; 
        m.min_val_x3 = val_x3_min;
        m.max_abs_val_diff_x3 = val_x3_diff_max;
        m.max_rel_val_diff_x3 = val_x3_diff_max / (val_x3_max - val_x3_min);

      }

      return m;
    }
  };

class qoi_x3_plugin : public libpressio_metrics_plugin {

  public:
    int begin_compress_impl(const struct pressio_data * input, struct pressio_data const * ) override {
      if(!input || !input->has_data()) return 0;
      input_data = pressio_data::clone(domain_manager().make_readable(domain_plugins().build("malloc"), *input));
      return 0;
    }
    int end_decompress_impl(struct pressio_data const*, struct pressio_data const* output, int ) override {
      if(!output || !output->has_data() || !input_data.has_data()) return 0;
      err_metrics = pressio_data_for_each<qoi_x3::qoi_x3_metrics>(input_data, domain_manager().make_readable(domain_plugins().build("malloc"), *output), qoi_x3::compute_metrics{});
      return 0;
    }

  struct pressio_options get_configuration_impl() const override {
    pressio_options opts;
    set(opts, "pressio:stability", "experimental");
    set(opts, "pressio:thread_safe", pressio_thread_safety_multiple);
    set(opts, "predictors:requires_decompress", true);
    set(opts, "predictors:invalidate", std::vector<std::string>{"predictors:error_dependent"});
    return opts;
  }


    struct pressio_options get_documentation_impl() const override {
      pressio_options opt;
    set(opt, "pressio:description", "qoi_x3");
    set(opt, "qoi_x3:max_val", "the maximum value squared");
    set(opt, "qoi_x3:min_val", "the minimum value squared");
    set(opt, "qoi_x3:max_abs_val_diff", "the maximum absolute difference of the squared values");
    set(opt, "qoi_x3:max_abs_rel_val_diff", "the maximum absolute relative difference of the squared values");
      return opt;
    }
    pressio_options get_metrics_results(pressio_options const &)  override {
      pressio_options opt;
      if(err_metrics) {
        set(opt, "qoi_x3:abx_error", (*err_metrics).max_abs_val_diff_x3);
        set(opt, "qoi_x3:rel_error", (*err_metrics).max_rel_val_diff_x3);
        set(opt, "qoi_x3:max_val", (*err_metrics).max_val_x3);
        set(opt, "qoi_x3:min_val", (*err_metrics).min_val_x3);
      } else {
        set_type(opt, "qoi_x3:abx_error", pressio_option_double_type);
        set_type(opt, "qoi_x3:rel_error", pressio_option_double_type);
        set_type(opt, "qoi_x3:max_val", pressio_option_double_type);
        set_type(opt, "qoi_x3:min_val", pressio_option_double_type);
      }
      return opt;
    }
    std::unique_ptr<libpressio_metrics_plugin> clone() override {
      return compat::make_unique<qoi_x3_plugin>(*this);
    }

  const char* prefix() const override {
    return "qoi_x3";
  }


  private:
  pressio_data input_data = pressio_data::empty(pressio_byte_dtype, {});
  compat::optional<qoi_x3::qoi_x3_metrics> err_metrics;

};

static pressio_register metrics_qoi_x3_plugin(metrics_plugins(), "qoi_x3", [](){ return compat::make_unique<qoi_x3_plugin>(); });
}}

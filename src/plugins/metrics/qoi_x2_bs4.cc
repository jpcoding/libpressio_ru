#include <algorithm>
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
#include <cstddef>
#include <memory>
#include <vector>
#include <iostream>  

namespace libpressio {
namespace qoi_x2_bs4 {
  struct qoi_x2_bs4_metrics {
    double abs_error;
    double rel_error;
    int block_size;
    double qoi_range; 
  };

  inline double qoi_fn(double x) {
    return x*x;
  }


template <class T>
void compute_qoi(T* data, T* dec_data, size_t num_elements){
  for(size_t i=0; i<num_elements; i++){
    data[i] = qoi_fn(data[i]);
    dec_data[i] = qoi_fn(dec_data[i]);
  }
}



    template <class T>
    std::vector<T> compute_average(T const * data, uint32_t n1, uint32_t n2, uint32_t n3, int block_size){
        uint32_t dim0_offset = n2 * n3;
        uint32_t dim1_offset = n3;
        uint32_t num_block_1 = (n1 - 1) / block_size + 1;
        uint32_t num_block_2 = (n2 - 1) / block_size + 1;
        uint32_t num_block_3 = (n3 - 1) / block_size + 1;
        std::vector<T> aggregated = std::vector<T>();
        uint32_t index = 0;
        T const * data_x_pos = data;
        for(int i=0; i<num_block_1; i++){
            int size_1 = (i == num_block_1 - 1) ? n1 - i * block_size : block_size;
            T const * data_y_pos = data_x_pos;
            for(int j=0; j<num_block_2; j++){
                int size_2 = (j == num_block_2 - 1) ? n2 - j * block_size : block_size;
                T const * data_z_pos = data_y_pos;
                for(int k=0; k<num_block_3; k++){
                    int size_3 = (k == num_block_3 - 1) ? n3 - k * block_size : block_size;
                    T const * cur_data_pos = data_z_pos;
                    int n_block_elements = size_1 * size_2 * size_3;
                    double sum = 0;
                    for(int ii=0; ii<size_1; ii++){
                        for(int jj=0; jj<size_2; jj++){
                            for(int kk=0; kk<size_3; kk++){
                                sum += *cur_data_pos;
                                cur_data_pos ++;
                            }
                            cur_data_pos += dim1_offset - size_3;
                        }
                        cur_data_pos += dim0_offset - size_2 * dim1_offset;
                    }
                    aggregated.push_back(sum / n_block_elements);
                    data_z_pos += size_3;
                }
                data_y_pos += dim1_offset * size_2;
            }
            data_x_pos += dim0_offset * size_1;
        }    
        return aggregated;
    }

    template<class T>
    T evaluate_L_inf(T const * data, T const * dec_data, uint32_t num_elements, bool normalized=true, bool verbose=false){
        T L_inf_error = 0;
        T L_inf_data = 0;
        for(int i=0; i<num_elements; i++){
            if(L_inf_data < fabs(data[i])) L_inf_data = fabs(data[i]);
            T error = data[i] - dec_data[i];
            if(L_inf_error < fabs(error)) L_inf_error = fabs(error);
        }
        return normalized ? L_inf_error / L_inf_data : L_inf_error;
    }




class qoi_x2_bs4_plugin : public libpressio_metrics_plugin {

  public:
    int begin_compress_impl(const struct pressio_data * input, struct pressio_data const * ) override {
      if(!input || !input->has_data()) return 0;
      input_data = pressio_data::clone(domain_manager().make_readable(domain_plugins().build("malloc"), *input));
      return 0;
    }
    int end_decompress_impl(struct pressio_data const*, struct pressio_data const* output, int ) override {
      if(!output || !output->has_data() || !input_data.has_data()) return 0;
      auto output_host = domain_manager().make_readable(domain_plugins().build("malloc"), *output);
      auto norm_dims = output_host.normalized_dims(4); 
      // compute_qoi<pressio_data_dtype(output)>(input_data.data(), output_host.data(), pressio_data_num_elements(output)); 
      // std::cout <<"num elements" <<pressio_data_num_elements(output)<< std::endl;
      // std::cout <<"dtype" <<output_host.dtype()<< std::endl; 
      size_t num_elements = pressio_data_num_elements(output); 

      int block_size = 4;


      if (output_host.dtype() ==pressio_float_dtype){
      using real= float;
      compute_qoi<real>((real*) input_data.data(),(real*) output_host.data(),  num_elements); 

      // std::cout << "dims" << norm_dims[2] << norm_dims[1] <<  norm_dims[0] << std::endl; 
      // std::cout << "max ele " << *std::max_element((real*) input_data.data(), (real*) input_data.data() + num_elements) << std::endl;
      if(norm_dims[2] ==0) norm_dims[2]= 1;
      auto avg_orig = qoi_x2_bs4::compute_average<real>((real*)(input_data.data()), norm_dims[2], norm_dims[1], norm_dims[0], block_size); 
      auto avg_decompressed = qoi_x2_bs4::compute_average<real>((real*) (output_host.data()), norm_dims[2], norm_dims[1], norm_dims[0], block_size); 
      auto L_inf_error = qoi_x2_bs4::evaluate_L_inf(avg_orig.data(), avg_decompressed.data(), avg_orig.size(), 0, false); 
      auto minmax = std::minmax_element(avg_orig.begin(),avg_orig.end());
      // std::cout << "min: " << *minmax.first << " max: " << *minmax.second << std::endl;
      // std::cout << "L_inf_error: " << L_inf_error << std::endl;
      // std::cout << "rel error " << L_inf_error / (*minmax.second - *minmax.first) << std::endl;  
      auto qoi_range = *minmax.second - *minmax.first;
      if (qoi_range ==0)
      {
        qoi_range = 1; 
      }
      err_metrics.emplace(qoi_x2_bs4_metrics{(double) L_inf_error, (double)  L_inf_error/qoi_range, block_size, (double) qoi_range});
    }
    else if (output_host.dtype() ==pressio_double_dtype){
      if(norm_dims[2] ==0) norm_dims[2]= 1;
      compute_qoi<double>((double*) input_data.data(),(double*) output_host.data(),  num_elements);
      auto avg_orig = qoi_x2_bs4::compute_average<double>((double*)(input_data.data()), norm_dims[2], norm_dims[1], norm_dims[0], block_size);
      auto avg_decompressed = qoi_x2_bs4::compute_average<double>((double*) (output_host.data()), norm_dims[2], norm_dims[1], norm_dims[0], block_size);
      auto L_inf_error = qoi_x2_bs4::evaluate_L_inf(avg_orig.data(), avg_decompressed.data(), avg_orig.size(), 0, false);
      auto minmax = std::minmax_element(avg_orig.begin(),avg_orig.end());
      // std::cout << "min: " << *minmax.first << " max: " << *minmax.second << std::endl;
      // std::cout << "L_inf_error: " << L_inf_error << std::endl;
      // std::cout << "rel error" << L_inf_error / (*minmax.second - *minmax.first) << std::endl;
      auto qoi_range = *minmax.second - *minmax.first;
      if (qoi_range ==0)
      {
        qoi_range = 1;
      }
      err_metrics.emplace(qoi_x2_bs4_metrics{L_inf_error, L_inf_error/qoi_range, block_size, qoi_range});
    }
    else{
      return 0;
    }
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
    set(opt, "pressio:description", "qoi_x2_bs4");
    set(opt, "qoi_x2_bs4:max_val", "the maximum value squared");
    set(opt, "qoi_x2_bs4:min_val", "the minimum value squared");
    set(opt, "qoi_x2_bs4:max_abs_val_diff", "the maximum absolute difference of the squared values");
    set(opt, "qoi_x2_bs4:max_abs_rel_val_diff", "the maximum absolute relative difference of the squared values");
      return opt;
    }
    pressio_options get_metrics_results(pressio_options const &)  override {
      pressio_options opt;
      if(err_metrics) {
        set(opt, "qoi_x2_bs4:abx_error", (*err_metrics).abs_error);
        set(opt, "qoi_x2_bs4:rel_error", (*err_metrics).rel_error);
        set(opt, "qoi_x2_bs4:range", (*err_metrics).qoi_range);
        set(opt, "qoi_x2_bs4:blocksize", (*err_metrics).block_size);
      } else {
        set_type(opt, "qoi_x2_bs4:abx_error", pressio_option_double_type);
        set_type(opt, "qoi_x2_bs4:rel_error", pressio_option_double_type);
        set_type(opt, "qoi_x2_bs4:range", pressio_option_double_type);
        set_type(opt, "qoi_x2_bs4:blocksize", pressio_option_int32_type);
      }
      return opt;
    }
    std::unique_ptr<libpressio_metrics_plugin> clone() override {
      return compat::make_unique<qoi_x2_bs4_plugin>(*this);
    }

  const char* prefix() const override {
    return "qoi_x2_bs4";
  }


  private:
  pressio_data input_data = pressio_data::empty(pressio_byte_dtype, {});
  compat::optional<qoi_x2_bs4::qoi_x2_bs4_metrics> err_metrics;

};

static pressio_register metrics_qoi_x2_bs4_plugin(metrics_plugins(), "qoi_x2_bs4", [](){ return compat::make_unique<qoi_x2_bs4_plugin>(); });
}}

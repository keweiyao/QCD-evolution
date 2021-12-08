#ifndef INTERP_ND_H
#define INTERP_ND_H

#include <vector>
#include <cmath>

typedef double (*InterpFtype)(std::vector<double>);

// A light-weighted N-dimensional table class that supports
// 1) Setting and retriving data using linear or N-dim indices
// 2) Interpolating within the range of the grid
// 3) Only equally-spaced grid.
class table_nd{
private:
    // the dimensional D of the table and the size of a unit cube 2^D
    const int dim_, power_dim_;
    // total number of data points
    int total_size_;
    // Linear data arrat, grid step, grid minimum and maximum
    std::vector<double> data_, step_, grid_min_, grid_max_;
    // Shape of each dimension [N1, ... ND], total_size = N1*N2...*ND
    // and the conjugate block size [N2*...*ND, N3*...ND, ..., ND, 1]
    std::vector<int> shape_, conj_size_;
    InterpFtype interpF_;
    
public:
    // initializer, given dimension, shape, grid minimum and maximum
    table_nd():
    dim_(0), power_dim_(1), 
    shape_({}), grid_min_({}), grid_max_({}),
    interpF_(NULL)
    {};
    table_nd(int dim, std::vector<int> shape, std::vector<double> gridmin, std::vector<double> gridmax, InterpFtype interpF):
    dim_(dim), power_dim_(std::pow(2,dim)), 
    shape_(shape), grid_min_(gridmin), grid_max_(gridmax),
    interpF_(interpF)
    {
        step_.resize(dim_);
        for (int i=0; i<dim_; i++) step_[i] = (grid_max_[i]-grid_min_[i])/(shape_[i]-1);
        conj_size_.resize(dim_);
        conj_size_[dim_-1] = 1;
        for (int i=1; i<dim_; i++) conj_size_[dim_-i-1] = conj_size_[dim_-i] * shape[dim_-i];
        total_size_ = conj_size_[0] * shape[0];
        data_.resize(total_size_);
    }
    // Convert a N-dim index to a linear index
    int ArrayIndex2LinearIndex(std::vector<int> & index){
        int k = 0;
        for (int i=0; i<dim_; i++) k += index[i] * conj_size_[i];
        return k;
    }
    // Convert a linear index to a N-dim index
    std::vector<int> LinearIndex2ArrayIndex(int k){
        std::vector<int> index;
        int K = k;
        for (int i=0; i<dim_; i++) {
            index.push_back(int(floor(K/conj_size_[i])));
            K = K%conj_size_[i];
        }
        return index;
    }
    // Convert a N-dim index to the physical values of each independent variable
    std::vector<double> ArrayIndex2Xvalues(std::vector<int> index){
        std::vector<double> xvals;
        for (int i=0; i<dim_; i++) xvals.push_back(grid_min_[i]+step_[i]*index[i]);
        return xvals;
    }
    int dim() const {return dim_;}; // return dimension
    int size() const {return total_size_;}; // return total size
    std::vector<int> shape() const { return shape_; }    // return shape
    std::vector<double> step() const{ return step_; }    // return step
    // set an element with a N-dim index
    void set_with_index(std::vector<int> index, double value) { data_[ArrayIndex2LinearIndex(index)] = value; }
    // set an element with a linear index
    void set_with_linear_index(int k, double value) { data_[k] = value; }
    // get an element with a N-dim index
    double get_with_index(std::vector<int> index) { return data_[ArrayIndex2LinearIndex(index)]; }
    // get an element with a linear index
    double get_with_linear_index(int k) { return data_[k]; }
    // Interpolating at a given array of physical varaibles X = [x1, x2, ..., xD]
    double interpolate(std::vector<double> Xinput) {
       // First, find the unit cube that contains the input point
       std::vector<int> start_index; start_index.resize(dim_);
       std::vector<double> w; w.resize(dim_);
       for(auto i=0; i<dim_; i++) {
           double x = (Xinput[i] - grid_min_[i])/step_[i];
           x = std::min(std::max(x, 0.), shape_[i]-1.);
           size_t nx = size_t(std::floor(x));
           double rx = x - nx;
           start_index[i] = nx;
           w[i] = rx;
       }
       // Second, use linear interpolation within the unit cube
       // interpolation is done over Grid/interpF and then multiple interpF
       // at the interped input
       std::vector<int> index(dim_);
       std::vector<double> indexInput(dim_);
       double result = 0.0;
       for (auto i=0; i<power_dim_; i++) {
           auto W = 1.0;
           for (auto j=0; j<dim_; j++) {
               index[j] = start_index[j] + ((i & ( 1 << j )) >> j);
               W *= (index[j]==start_index[j])?(1.-w[j]):w[j];
           }

           for (auto j=0; j<dim_; j++) 
               indexInput[j] = grid_min_[j] + index[j]*step_[j];

           result += data_[ArrayIndex2LinearIndex(index)]
                   / interpF_(indexInput) * W;
       }
       return result*interpF_(Xinput);
    }
};

#endif

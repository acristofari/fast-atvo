#ifndef __FAST_ATVO_H_INCLUDED__
#define __FAST_ATVO_H_INCLUDED__
#include <vector>
#include <fstream>

class Graph;

struct fast_atvo_options {
    unsigned int ws_size,out_it;
    unsigned short int verbosity;
    double p_exp,lb,ub,perc_at_bounds;
};

class Fast_atvo {
private:
    const Graph *gr;

    // parameters for termination of local minimization (to be set in the constructor)
    unsigned int max_it,max_n_f,max_n_g;
    double eps_opt,min_gd,min_norm_proj_d,min_stepsize,min_decrease_f;

    unsigned int n,m,z,z_nm,it,it_bh,k,it_nm,n_ws,n_new_act_l,n_new_act_u;
    unsigned int n_f,n_g,seed,n_ws_max,n_ws_max_k;
    unsigned short int flag,verbosity;
    double p_exp,modularity,l,u,f,f_best_local,f_prev,f_w,fv,stepsize,sq_norm_g_ws;
    double norm_proj_d,sup_norm_proj_g,gd,f_first,delta_f0;
    bool ls,f_computed,gd_exit,dir_exit,stepsize_exit,f_exit,is_first_linesearch,is_ws_full;
    std::vector<std::vector<double>> gr_inst_val_g;
    std::vector<unsigned int> ind_ws;
    std::vector<unsigned short int> c;
    std::vector<double> x,x_best_local,x_original,g,g_best_local,g_old,v,w,d;
    std::vector<bool> is_in_ws;
    std::ofstream file_output;
    
    const short unsigned int solve_locally();
    void grad(),grad_fast(),swap_vars(const std::vector<double>&),main_prints();
    void compute_working_set(),update_w(),restart(),linesearch(bool);
    void compute_communities(const std::vector<double>&);
    inline void compute_sup_norm_proj_g();
    bool converged(),check_stop();
    inline double sign(const double&);

public:
    Fast_atvo(const Graph*, const std::vector<double>&, const fast_atvo_options&);
    void solve();
    const std::vector<unsigned short int>& get_communities();
    double get_modularity();
    const std::vector<double>& get_x();
};

#endif
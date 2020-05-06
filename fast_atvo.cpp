// -------------------------------------------------------------------------
//
// This file is part of FAST-ATVO, which is a software for community
// detection in an undirected graph with non-negative weights.
// See 'README.txt' to see how to use FAST-ATVO.
//
// -------------------------------------------------------------------------
//
// Reference paper:
// A. Cristofari, F. Rinaldi, F. Tudisco (2020). Total variation based
// community detection using a nonlinear optimization approach. SIAM Journal
// on Applied Mathematics, to appear
//
// -------------------------------------------------------------------------
//
// Authors:
// Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
// Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
// Francesco Tudisco (e-mail: francesco.tudisco@gssi.it)
//
// Last update:
// May 6th, 2020
//
// Copyright 2019-2020 Andrea Cristofari, Francesco Rinaldi, Francesco
// Tudisco.
//
// Licensing:
// This file is part of FAST-ATVO.
// FAST-ATVO is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// FAST-ATVO is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with FAST-ATVO. If not, see <http://www.gnu.org/licenses/>.
//
// -------------------------------------------------------------------------

#include <numeric>
#include <algorithm>
#include <iostream>
#include <random>
#include <math.h>
#include "fast_atvo.h"
#include "graph.h"


// constructor
//-------------------------------------------------------------------------------------
Fast_atvo::Fast_atvo(const Graph *p, const std::vector<double>& x0, const alg_options& alg_opts) {
            
    // parameters for termination of local minimization
    // (see description of 'solve_locally' function below)
    //----------------------------------------------------------------------
    eps_opt = 1e-1;
    min_gd = 1e-6;
    min_norm_proj_d = 1e-7;
    min_stepsize = 1e-20;
    min_decrease_f = 1e-2;
    max_it = 1000;
    max_n_f = 1000;
    max_n_g = 1000;
    //----------------------------------------------------------------------

    // non-monotone parameters for local minimization
    //----------------------------------------------------------------------
    m = 100;
    z = 20;
    //----------------------------------------------------------------------

    seed = 1;
    srand(seed);
    gr = p;
    n = gr->n;
    p_exp = alg_opts.p_exp;
    n_ws_max = (unsigned int) alg_opts.ws_size;
    it_bh = (unsigned int) alg_opts.out_it;
    l = alg_opts.lb;
    u = alg_opts.ub;
    verb = alg_opts.verbosity;
    if (u-l < 2e0) {
        eps_opt *= (u-l)/2e0;
        min_norm_proj_d *= (u-l)/2e0;
        min_decrease_f *= (u-l)/2e0;
    }
    x.resize(n);
    g.resize(n);
    if (alg_opts.perc_at_bounds >= 1e0) {
        for (unsigned int i=0; i<n; i++) {
            x[i] = x0[gr->ind_n_to_n_original[i]] < 0e0 ?
                   l : (x0[gr->ind_n_to_n_original[i]] > 0e0 ? u : 0e0);
        }
    } else if (alg_opts.perc_at_bounds > 0e0) {
        unsigned int count_l,count_u;
        std::vector<unsigned int> ind_vars;
        ind_vars.resize(n);
        count_l = 0;
        count_u = n - 1;
        for (unsigned int i=0; i<n; i++) {
            x[i] = std::max(l,std::min(x0[gr->ind_n_to_n_original[i]],u));
            if (x0[i] < 0e0) {
                ind_vars[count_l] = i;
                count_l++;
            } else if (x0[i] > 0e0) {
                ind_vars[count_u] = i;
                count_u--;
            }
        }
        count_u = n + 1 - count_u;
        shuffle(ind_vars.begin(),ind_vars.begin()+count_l,std::default_random_engine(seed++));
        shuffle(ind_vars.end()-count_u,ind_vars.end(),std::default_random_engine(seed++));
        for (unsigned int i=0; i<std::min((unsigned int)ceil(alg_opts.perc_at_bounds*count_l),(unsigned int)count_l); i++) {
            x[i] = l;
        }
        for (unsigned int i=n-std::min((unsigned int)ceil(alg_opts.perc_at_bounds*count_u),(unsigned int)count_u); i<n; i++) {
            x[i] = u;
        }
    } else {
        for (unsigned int i=0; i<n; i++) {
            x[i] = std::max(l,std::min(x0[gr->ind_n_to_n_original[i]],u));
        }
    }
    
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Fast_atvo::solve() {
    if (verb >= 1) {
        std::cout.precision(4);
        std::cout.setf(std::ios::scientific,std::ios::floatfield);
        file_output.open("iteration_history.txt",std::ios::trunc);
        file_output.precision(4);
        file_output.setf(std::ios::scientific,std::ios::floatfield);
        std::cout << "number of variables = " << n << " ("
                  << gr->n_original-n << " were removed in pre-process)\n\n\n*** outer iteration 1 ***";
        file_output << "number of variables = " << n << " ("
                    << gr->n_original-n << " were removed in pre-process)\n\n\n*** outer iteration 1 ***";
    }
    gr_inst_val_g = gr->inst.val;
    for (unsigned int i=0; i<gr->inst.n_rows; i++) {
        for (unsigned int j=0; j<gr->inst.val[i].size(); j++) {
            gr_inst_val_g[i][j] -= gr->degree[gr->inst.node_index[i]]*gr->degree[gr->inst.ind_row_col[i][j]]/gr->volume;
        }
    }
    solve_locally();
    compute_communities(x);
    if (it_bh > 1) {
        std::vector<double> x_best = x;
        std::vector<unsigned short int> c_best = c;
        double modularity_best = modularity;
        double eps_opt_orig = eps_opt;
        eps_opt = std::min(9e-1*(u-l),pow(15e-1,floor((double)(it_bh)/2e0))*eps_opt_orig);
        for (unsigned int i=1; i<it_bh; i++) {
            swap_vars(x_best);
            if (verb >= 1) {
                std::cout << "\n\n\n*** outer iteration " <<  i + 1 << " ***";
                file_output << "\n\n\n*** outer iteration " <<  i + 1 << " ***";
            }
            eps_opt = std::max(eps_opt_orig,eps_opt/15e-1);
            solve_locally();
            compute_communities(x);
            if (modularity > modularity_best) {
                x_best = x;
                c_best = c;
                modularity_best = modularity;
            }
        }
        x = x_best;
        c = c_best;
        modularity = modularity_best;
    }
    if (verb >= 1) {
        std::cout << "\n\nWARNING: using 'verb = 0' may be faster\n\n";
        file_output << "\n\nWARNING: using 'verb = 0' may be faster\n";
        file_output.close();
    }
    // clear vectors to free memory (only 'c' and 'x' are maintained)
    gr_inst_val_g.clear();
    gr_inst_val_g.shrink_to_fit();
    ind_ws.clear();
    ind_ws.shrink_to_fit();
    x_best_local.clear();
    x_best_local.shrink_to_fit();
    g.clear();
    g.shrink_to_fit();
    g_best_local.clear();
    g_best_local.shrink_to_fit();
    g_old.clear();
    g_old.shrink_to_fit();
    v.clear();
    v.shrink_to_fit();
    w.clear();
    w.shrink_to_fit();
    d.clear();
    d.shrink_to_fit();
    is_in_ws.clear();
    is_in_ws.shrink_to_fit();
    if (!x_original.empty()) {
        x_original.clear();
        x_original.shrink_to_fit();
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Fast_atvo::swap_vars(const std::vector<double>& x_best) {

    unsigned int count_l,count_u;
    std::vector<unsigned int> ind_vars;
    std::minstd_rand0 eng(seed);
    const double double_rand_max = (double) RAND_MAX;
    const double r = 75e-2; // percentage of variables to be swapped

    ind_vars.resize(n);
    count_l = 0;
    count_u = n - 1;

    for (unsigned int i=0; i<n; i++) {
        if (x_best[i] < 0e0) {
            ind_vars[count_l++] = i;
        } else if (x_best[i] > 0e0) {
            ind_vars[count_u--] = i;
        } else if (double(rand())/double_rand_max < 5e-1) {
            ind_vars[count_l++] = i;
        } else {
            ind_vars[count_u--] = i;
        }
    }

    count_u = n - count_u + 1;
    shuffle(ind_vars.begin(),ind_vars.begin()+count_l,std::default_random_engine(seed=eng()));
    shuffle(ind_vars.end()-count_u,ind_vars.end(),std::default_random_engine(seed=eng()));
    x = x_best;
    for (unsigned int i=0; i<std::min((unsigned int)ceil(r*count_l),(unsigned int)count_l); i++) {
        x[i] = u;
    }
    for (unsigned int i=n-std::min((unsigned int)ceil(r*count_u),(unsigned int)count_u); i<n; i++) {
        x[i] = l;
    }

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const short unsigned int Fast_atvo::solve_locally() {

    // Optimization solver
    //
    // The output is an integer describing the exit condition:
    //   0 if the sup-norm of x - p[x-g(x)] <= eps_opt*min(2,(u-l))/2, where
    //     'x' is the solution found by the algorithm, 'l' and 'u' are the lower
    //     and the upper bound on the variables, respectively, 'p[.]' is the
    //     projection operator onto the feasible set and 'g(.)' is the gradient
    //     of the objective function
    //     (default value of 'eps_opt' = 1e-1),
    //   1 if the directional derivative (in absolute value) of the objective
    //     function along the search direction <= min_gd
    //     (default value of 'min_gd' = 1e-15),
    //   2 if the norm of the projected search direction <=
    //     min_norm_proj_d*min(2,(u-l))/2, where 'l' and 'u' are the lower and
    //     the upper bound on the variables, respectively
    //     (default value of 'min_norm_proj_d' = 1e-7),
    //   3 if the stepsize <= min_stepsize
    //     (default value of 'min_stepsize' = 1e-20),
    //   4 if the the objective decrease (in percentage) <=
    //     min_decrease_f*min(2,(u-l))/2, where 'l' and 'u' are the lower and
    //     the upper bound on the variables, respectively
    //     (default value of 'min_decrease_f' = 1e-2),
    //   5 if the number of iterations >= max_it
    //     (default value of 'max_it' = 1000),
    //   6 if the number of function evaluations >= max_n_f
    //     (default value of 'max_n_f' = 1000),
    //   7 if the number of gradient evaluations >= max_n_g
    //     (default value of 'max_n_g' = 1000)

    unsigned int j;
    double delta0_dir,delta_dir,beta_dir,tmp,c_bb,sy_bb;
    std::vector<double> s_bb,y_bb;
    bool warn_bb,is_restarted;
    const double c_bb_min = 1e-10;
    const double c_bb_max = 1e10;

    // initialize vectors
    ind_ws.assign(n,0); // the indices of the variables included in the working set
                        // are in the first 'n_ws' positions
    s_bb.assign(n,0e0);
    y_bb.assign(n,0e0);
    d.assign(n,0e0);

    // compute the gradient
    grad();
    n_g = 0;
    
    // compute the objective function
    f = std::inner_product(g.begin(),g.end(),x.begin(),0e0)/p_exp;
    n_f = 0;

    x_best_local = v = x;
    f_best_local = f_w = f;
    g_best_local = g_old = g;

    w.assign(m,f);

    gd = std::numeric_limits<double>::lowest();
    
    // compute the sup-norm of the projected gradient
    compute_sup_norm_proj_g();
    
    delta0_dir = delta_dir = 1e20;
    beta_dir = 99e-2; // reduction factor of 'delta_dir' (must be >=0 and <1)

    z_nm = std::min(n,z);

    // initialize counters
    it = k = it_nm = 0;

    n_ws_max_k = 2;
    
    gd_exit = false;
    dir_exit = false;
    stepsize_exit = false;
    f_exit = false;
    f_computed = true;
    is_first_linesearch = true;

    //-----------------
    // START MAIN LOOP
    //-----------------

    while (!converged()) {

        is_restarted = false;
        compute_working_set();

        if (verb > 1) {
            std::cout << "\n\n--- iteration details ---\n\nsize of the working set = " << n_ws;
            file_output << "\n\n--- iteration details ---\n\nsize of the working set = " << n_ws;
        }

        sq_norm_g_ws = 0e0;
        for (unsigned int i=0; i<n_ws; i++) {
            sq_norm_g_ws += g[ind_ws[i]]*g[ind_ws[i]];
        }

        // function control
        if (it_nm >= z_nm) {
            z_nm = std::min(z_nm+n,z);
            if (!f_computed) {
                f = std::inner_product(g.begin(),g.end(),x.begin(),0e0)/p_exp; // evaluate f(x)
                n_f++;
            }
            if (f >= f_w) {
                if (verb > 1) {
                    std::cout << "\nfunction control not satisfied (f = " << f << ")\n\nrestart from the best point";
                    file_output << "\nfunction control not satisfied (f = " << f << ")\n\nrestart from the best point";
                }
                restart();
            } else {
                if (verb > 1) {
                    std::cout << "\nfunction control satisfied (f = " << f << ")";
                    file_output << "\nfunction control satisfied (f = " << f << ")";
                }
                update_w();
                f_computed = true;
            }
        }

        if (!is_restarted) {

            // compute the search direction
            // N.B. 'v' is the previous iterate
            if (k >= 2) {
                sy_bb = 0e0;
                for (unsigned int i=0; i<n_ws; i++) {
                    j = ind_ws[i];
                    s_bb[i] = x[j] - v[j];
                    y_bb[i] = g[j] - g_old[j];
                    sy_bb += s_bb[i]*y_bb[i];
                }
                if (sy_bb > 0e0) {
                    c_bb = std::inner_product(s_bb.begin(),s_bb.begin()+n_ws,s_bb.begin(),0e0)/sy_bb;
                    if (c_bb > c_bb_min) {
                        c_bb = std::min(c_bb,c_bb_max);
                        warn_bb = false;
                    } else {
                        c_bb = sy_bb/std::inner_product(y_bb.begin(),y_bb.begin()+n_ws,y_bb.begin(),0e0);
                        if (c_bb > c_bb_min) {
                            c_bb = std::min(c_bb,c_bb_max);
                            warn_bb = false;
                        } else {
                            c_bb = c_bb_min;
                            warn_bb = true;
                        }
                    }
                } else {
                    tmp = 0e0; // ||x_{W^k}||^2
                    for (unsigned int i=0; i<n_ws; i++) {
                        tmp += x[ind_ws[i]]*x[ind_ws[i]];
                    }
                    c_bb =  std::min(c_bb_max,std::max(1e0,sqrt(tmp/sq_norm_g_ws)));
                    warn_bb = true;
                }
            } else {
                tmp = 0e0; // ||x_{W^k}||^2
                for (unsigned int i=0; i<n_ws; i++) {
                    tmp += x[ind_ws[i]]*x[ind_ws[i]];
                }
                c_bb =  std::min(c_bb_max,std::max(1e0,sqrt(tmp/sq_norm_g_ws)));
                warn_bb = true;
            }
            for (unsigned int i=0; i<n_ws; i++) {
                d[ind_ws[i]] = -c_bb*g[ind_ws[i]];
            }
            gd = -c_bb*sq_norm_g_ws;

            if (verb > 1) {
                std::cout << "\ndirectional derivative = " << gd;
                file_output << "\ndirectional derivative = " << gd;
            }

            // check if the directional derivative is sufficiently negative
            if (gd < -min_gd) {

                // check if line search must be performed
                if (is_first_linesearch || warn_bb) {
                    ls = true;
                } else {
                    ls = false;
                }

                v = x;

                // try accepting the unit stepsize or prepare for the line search
                if (!ls) {

                    // set v = x and x = p(x + d)
                    norm_proj_d = 0e0;
                    for (unsigned int i=0; i<n_ws; i++) {
                        j = ind_ws[i];
                        x[j] = std::max(l,std::min(x[j]+d[j],u));
                        norm_proj_d += (x[j]-v[j])*(x[j]-v[j]);
                    }
                    norm_proj_d = sqrt(norm_proj_d);

                    if (verb > 1) {
                        std::cout << "\nnorm of the projected direction = " << norm_proj_d;
                        file_output << "\nnorm of the projected direction = " << norm_proj_d;
                    }

                    // check if the projected direction if sufficiently large
                    if (norm_proj_d > min_norm_proj_d) {
                        if (norm_proj_d <= delta_dir) { // test satisfied -> unit stepsize accepted
                            delta_dir *= beta_dir;
                            g_old = g;
                            // compute the gradient and the sup-norm of the projected gradient
                            if (n_ws <= double(n-1)/3e0) {
                                grad_fast();
                            } else {
                                grad();
                            }
                            n_g++;
                            compute_sup_norm_proj_g();
                            f_computed = false;
                            it_nm++;
                            if (verb > 1) {
                                std::cout << "\nunit stepsize accepted without computing f";
                                file_output << "\nunit stepsize accepted without computing f";
                            }
                        } else if (!f_computed) { // check the objective function
                            // restore x
                            for (unsigned int i=0; i<n_ws; i++) {
                                j = ind_ws[i];
                                tmp = v[j];
                                v[j] = x[j];
                                x[j] = tmp;
                            }
                            f = std::inner_product(g.begin(),g.end(),x.begin(),0e0)/p_exp; // evaluate f(x)
                            n_f++;
                            if (f < f_w) { // objective function decreased -> line search
                                update_w();
                                f_computed = true;
                                ls = true;
                                // set again v = x and x = p(x + d)
                                for (unsigned int i=0; i<n_ws; i++) {
                                    j = ind_ws[i];
                                    tmp = v[j];
                                    v[j] = x[j];
                                    x[j] = tmp;
                                }
                                if (verb > 1) {
                                    std::cout << "\nfunction control satisfied (f = " << f << ")";
                                    file_output << "\nfunction control satisfied (f = " << f << ")";
                                }
                            } else { // objective function not decreased -> restart
                                if (verb > 1) {
                                    std::cout << "\nfunction control not satisfied (f = " << f
                                              << ")\n\nrestart from the best point";
                                    file_output << "\nfunction control not satisfied (f = " << f
                                                << ")\n\nrestart from the best point";
                                }
                                restart();
                            }
                        } else {
                            ls = true;
                        }
                    } else {
                        // restore x
                        for (unsigned int i=0; i<n_ws; i++) {
                            x[ind_ws[i]] = v[ind_ws[i]];
                        }
                        dir_exit = true;
                    }

                } else if (!f_computed) {

                    f = std::inner_product(g.begin(),g.end(),x.begin(),0e0)/p_exp;  // evaluate f(x)
                    n_f++;
                    if (f < f_w) { // objective function decreased -> line search
                        // set v = x and x = p(x + d)
                        f_computed = true;
                        norm_proj_d = 0e0;
                        for (unsigned int i=0; i<n_ws; i++) {
                            j = ind_ws[i];
                            x[j] = std::max(l,std::min(x[j]+d[j],u));
                            norm_proj_d += (x[j]-v[j])*(x[j]-v[j]);
                        }
                        norm_proj_d = sqrt(norm_proj_d);
                        if (verb > 1) {
                            std::cout << "\nnorm of the projected direction = " << norm_proj_d;
                            file_output << "\nnorm of the projected direction = " << norm_proj_d;
                        }
                        if (norm_proj_d > min_norm_proj_d) {
                            update_w();
                            if (verb > 1) {
                                std::cout << "\nfunction control satisfied (f = " << f << ")";
                                file_output << "\nfunction control satisfied (f = " << f << ")";
                            }
                        } else {
                            // restore x
                            for (unsigned int i=0; i<n_ws; i++) {
                                x[ind_ws[i]] = v[ind_ws[i]];
                            }
                            dir_exit = true;
                            ls = false; // to skip the line search and exit the while loop
                        }
                    } else { // objective function not decreased -> restart
                             // (first check if the projected direction is sufficiently large)
                        norm_proj_d = 0e0;
                        for (unsigned int i=0; i<n_ws; i++) {
                            j = ind_ws[i];
                            tmp = x[j] - std::max(1e0,std::min(x[j]+d[j],1e0));
                            norm_proj_d += tmp*tmp;
                        }
                        norm_proj_d = sqrt(norm_proj_d);
                        if (verb > 1) {
                            std::cout << "\nnorm of the projected direction = " << norm_proj_d;
                            file_output << "\nnorm of the projected direction = " << norm_proj_d;
                        }
                        if (norm_proj_d > min_norm_proj_d) {
                            if (verb > 1) {
                                std::cout << "\npoint not accepted (f = " << f
                                          << ")\n\nrestart from the best point";
                                file_output << "\npoint not accepted (f = " << f
                                            << ")\n\nrestart from the best point";
                            }
                            restart();
                        } else {
                            dir_exit = true;
                            f_computed = true;
                            ls = false; // to skip the line search and exit the while loop
                        }                            
                    }

                } else {

                    // set v = x and x = p(x + d)
                    norm_proj_d = 0e0;
                    for (unsigned int i=0; i<n_ws; i++) {
                        j = ind_ws[i];
                        x[j] = std::max(l,std::min(x[j]+d[j],u));
                        norm_proj_d += (x[j]-v[j])*(x[j]-v[j]);
                    }
                    norm_proj_d = sqrt(norm_proj_d);
                    if (verb > 1) {
                        std::cout << "\nnorm of the projected direction = " << norm_proj_d;
                        file_output << "\nnorm of the projected direction = " << norm_proj_d;
                    }
                    if (norm_proj_d <= min_norm_proj_d) {
                        // restore x
                        for (unsigned int i=0; i<n_ws; i++) {
                            x[ind_ws[i]] = v[ind_ws[i]];
                        }
                        dir_exit = true;
                        ls = false; // to skip the line search and exit the while loop
                    }

                }

                // line search
                if (ls && !is_restarted) {
                    if (verb >= 2) {
                        std::cout << "\nline search";
                        file_output << "\nline search";
                    }
                    if (is_first_linesearch) {
                        f_first = f;
                        delta_f0 = 0e0;
                    }
                    f_prev = f;
                    g_old = g;
                    linesearch(n_ws <= double(n-1)/3e0); // N.B. 'f_computed' = true at this point
                    if (!stepsize_exit) {
                        n_g++;                        
                        compute_sup_norm_proj_g();
                        
                        if (verb > 1) {
                            std::cout << "\nstepsize = " << stepsize;
                            file_output << "\nstepsize = " << stepsize;
                        }
                        if ((f_prev-f)/(fabs(f_prev)+1e0)>min_decrease_f || !is_ws_full) {
                            if (is_first_linesearch) {
                                delta_f0 = f_first - f;
                                delta_dir = delta0_dir*stepsize*norm_proj_d;
                                is_first_linesearch = false;
                            }
                            update_w();
                        } else {
                            f_exit = true;
                        }
                    } else {
                        // restore x and f(x)
                        for (unsigned int i=0; i<n_ws; i++) {
                            x[ind_ws[i]] = v[ind_ws[i]];
                        }
                        f = fv;
                    }
                }

            } else {

                gd_exit = true;

            }

        }

        it++;
        k++;
        n_ws_max_k = std::min(2*n_ws_max_k,n_ws_max);

    }

    return flag;

}
//-------------------------------------------------------------------------------------




// other functions
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
bool Fast_atvo::converged() {

    if ( (sup_norm_proj_g > eps_opt) && (!gd_exit) && (!dir_exit) && (!stepsize_exit) &&
        (!f_exit) && (it<max_it) && (n_f<max_n_f) && (n_g<max_n_g) ) {
        if (verb > 0) {
            main_prints();
        }
        return false;
    } else {
        if (!f_computed) {
            f = std::inner_product(g.begin(),g.end(),x.begin(),0e0)/p_exp;
            n_f++;
            f_computed = true;
        }
        if (sup_norm_proj_g <= eps_opt) {
            if (verb > 0) {
                main_prints();
                std::cout << "\n\n================================================================================"
                          << "\noptimality condition satisfied: sup-norm of the projected gradient <= " << eps_opt
                          << "\n================================================================================";
                file_output << "\n\n================================================================================"
                            << "\noptimality condition satisfied: sup-norm of the projected gradient <= " << eps_opt
                            << "\n================================================================================";
            }
            if (f <= f_best_local) {
                flag = 0;
                return true;
            } else {
                return check_stop();
            }
        } else if (gd_exit) {
            if (verb > 0) {
                std::cout << "\n\ndirectional derivative not sufficiently negative";
                file_output << "\n\ndirectional derivative not sufficiently negative";
            }
            if (f <= f_best_local) {
                it--;
                flag = 1;
                if (verb > 1) {
                    std::cout << ", algorithm stopped";
                    file_output << ", algorithm stopped";
                }
                return true;
            } else {
                return check_stop();
            }
        } else if (dir_exit) {
            if (verb > 0) {
                std::cout << "\n\nprojected direction too small";
                file_output << "\n\nprojected direction too small";
            }
            if (f <= f_best_local) {
                it--;
                flag = 2;
                if (verb > 1) {
                    std::cout << ", algorithm stopped";
                    file_output << ", algorithm stopped";
                }
                return true;
            } else {
                return check_stop();
            }
        } else if (stepsize_exit) {
            if (verb > 0) {
                std::cout << "\n\nstepsize too small";
                file_output << "\n\nstepsize too small";
            }
            if (f <= f_best_local) {
                it--;
                flag = 3;
                if (verb > 1) {
                    std::cout << ", algorithm stopped";
                    file_output << ", algorithm stopped";
                }
                return true;
            } else {
                return check_stop();
            }
        } else if (f_exit) {
            if (verb > 0) {
                std::cout << "\n\nobjective decrease too small";
                file_output << "\n\nobjective decrease too small";
            }
            if (f <= f_best_local) {
                flag = 4;
                if (verb > 1) {
                    std::cout << ", algorithm stopped";
                    file_output << ", algorithm stopped";
                }
                return true;
            } else {
                return check_stop();
            }
        } else {
            if (f > f_best_local) {
                x = x_best_local;
            }
            if (it >= max_it) {
                flag = 5;
                if (verb > 0) {
                    main_prints();
                    std::cout << "\n\ntoo many iterations, algorithm stopped";
                    file_output << "\n\ntoo many iterations, algorithm stopped";
                }
            } else if (n_f >= max_n_f) {
                flag = 6;
                if (verb > 0) {
                    main_prints();
                    std::cout << "\n\ntoo many function evaluations, algorithm stopped";
                    file_output << "\n\ntoo many function evaluations, algorithm stopped";
                }
            } else if (n_g >= max_n_g) {
                flag = 7;
                if (verb > 0) {
                    main_prints();
                    std::cout << "\n\ntoo many gradient evaluations, algorithm stopped";
                    file_output << "\n\ntoo many gradient evaluations, algorithm stopped";
                }
            }
            return true;
        }
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
bool Fast_atvo::check_stop() { // it returns true if the algorithm must stop, false otherwise
    if ( (it<max_it) && (n_f<max_n_f) && (n_g<max_n_g) ) {
        restart();
        if (verb > 0) {
            if (verb > 1) {
                std::cout << "\n\nrestart from the best point";
                file_output << "\n\nrestart from the best point";
            }
            main_prints();
        }
        return false;
    } else {
        x = x_best_local;
        if (it >= max_it) {
            flag = 5;
            if (verb > 0) {
                std::cout << "\ntoo many iterations, algorithm stopped";
                file_output << "\ntoo many iterations, algorithm stopped";
            }
        } else if (n_f >= max_n_f) {
            flag = 6;
            if (verb > 0) {
                std::cout << "\ntoo many function evaluations, algorithm stopped";
                file_output << "\ntoo many function evaluations, algorithm stopped";
            }
        } else { // n_g >= max_n_g
            flag = 7;
            if (verb > 0) {
                std::cout << "\ntoo many gradient evaluations, algorithm stopped";
                file_output << "\ntoo many gradient evaluations, algorithm stopped";
            }
        }
        return true;
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Fast_atvo::main_prints() {
    std::cout << "\n\n--------------------------------------------------------\n\n"
              << "iteration " << it
              << "\n\nbest f = " << f_best_local
              << "\nsup-norm of the projected gradient at the current point = " << sup_norm_proj_g
              << "\nnumber of function evaluations = " << n_f
              << "\nnumber of gradient evaluations = " << n_g;
    file_output << "\n\n--------------------------------------------------------\n\n"
                << "iteration " << it
                << "\n\nbest f = " << f_best_local
                << "\nsup-norm of the projected gradient at the current point = " << sup_norm_proj_g
                << "\nnumber of function evaluations = " << n_f
                << "\nnumber of gradient evaluations = " << n_g;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
inline void Fast_atvo::compute_sup_norm_proj_g() { // compute the sup-norm of (x - p[x-g(x)])
    sup_norm_proj_g = 0e0;
    for (unsigned int i=0; i<n; i++) {
        sup_norm_proj_g = std::max(sup_norm_proj_g,fabs(x[i]-std::max(l,std::min(x[i]-g[i],u))));
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Fast_atvo::compute_working_set() {
    unsigned int j,h;
    std::minstd_rand0 eng(seed);
    
    // non-active set estimate
    // (variables with zero derivative are ignored, since the corresponding
    // components of the search direction would be zeros)
    n_ws = 0;
    is_in_ws.assign(n,false);
    for (unsigned int i=0; i<n; i++) {
        if ( !(x[i]<=l && g[i]>0e0) && !(x[i]>=u && g[i]<0e0) && g[i]!=0e0) {
            ind_ws[n_ws] = i;
            is_in_ws[i] = true;
            n_ws++;
        }
    }
     
    if (n_ws > n_ws_max_k) { // maximum size allowed for the working set exceeded
        is_in_ws.assign(n,false);
        shuffle(ind_ws.begin(),ind_ws.begin()+n_ws,std::default_random_engine(seed=eng()));
        double sup_norm_proj_g_ws = fabs(x[0]-std::max(l,std::min(x[0]-g[0],u)));
        for (unsigned int i=1; i<n_ws_max_k; i++) {
            h = ind_ws[i];
            is_in_ws[h] = true;
            sup_norm_proj_g_ws = std::max(sup_norm_proj_g_ws,fabs(x[h]-std::max(l,std::min(x[h]-g[h],u))));
        }
        if (sup_norm_proj_g_ws < sup_norm_proj_g) {
            unsigned int i_max = ind_ws[0];
            double tmp = sup_norm_proj_g_ws;
            j = n_ws_max_k;
            while (sup_norm_proj_g_ws<sup_norm_proj_g && j<n) {
                h = ind_ws[j];
                tmp = fabs(x[h]-std::max(l,std::min(x[h]-g[h],u)));
                if (tmp > sup_norm_proj_g_ws) {
                    i_max = h;
                    sup_norm_proj_g_ws = tmp;
                }
                j++;
            }
            ind_ws[0] = i_max;
            is_in_ws[i_max] = true;
        } else {
            is_in_ws[ind_ws[0]] = true;
        }
        n_ws = n_ws_max_k;
        is_ws_full = false;
    } else {
        is_ws_full = true;
    }

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Fast_atvo::update_w() { // update the vector of reference values
    if (f < f_best_local) {
        x_best_local = x;
        f_best_local = f;
        g_best_local = g;
    }
    f_w = f;
    for (unsigned int i=m-1; i>0; i--) {
        w[i] = w[i-1];
        f_w = std::max(f_w,w[i]);
    }
    w[0] = f;
    it_nm = 0;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Fast_atvo::restart() {
    x = v = x_best_local;
    f = f_w = f_best_local;
    g = g_old = g_best_local;
    m = (m+4)/5;
    z_nm = std::min(n,z);
    w.assign(m,f_w);
    k = it_nm = 0;
    n_ws_max_k = 2;
    //delta0_dir *= 1e-1;
    gd_exit = false;
    dir_exit = false;
    stepsize_exit = false;
    f_exit = false;
    f_computed = true;
    is_first_linesearch = true;
    compute_sup_norm_proj_g();
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Fast_atvo::linesearch(bool fast_f) {
    unsigned int j;
    double f_ref,aa,a1;
    bool red_param_computed;
    const double gamma = 1e-3;
    const double delta = 5e-1;
    const double f_dec = 1e6;

    // N.B. 'v' is the current point, 'x' is the point obtained by using
    // the search direction with unit stepsize (and projecting onto the box)
    // and 'f' is f(v)

    red_param_computed = false;

    if (k >= 2) {
        stepsize = 1e0;
    } else {
        if (k == 0) {
            stepsize = std::max(std::pow(delta,8e0),4e0*min_stepsize);
        }
        for (unsigned int i=0; i<n_ws; i++) {
            j = ind_ws[i];
            x[j] = std::max(l,std::min(v[j]+stepsize*d[j],u));
        }
    }

    fv = f;
    if (fast_f) {
        grad_fast();
    } else {
        grad();
    }
    f = std::inner_product(g.begin(),g.end(),x.begin(),0e0)/p_exp;
    n_f++;

    if (f-f_first >= f_dec*delta_f0) {
        f_ref = fv;
    } else {
        f_ref = f_w;
    }

    while (f > f_ref+gamma*stepsize*gd) {

        // update the stepsize
        if ((f-fv)/std::min(1e12,std::max(1e-12,-gamma*stepsize*gd)) < 1e10) {
            stepsize *= delta;
        } else {
            if (!red_param_computed) {
                a1 = 1e0*std::max(1e0,sqrt(std::inner_product(v.begin(),v.end(),v.begin(),0e0)))/std::max(1e-15,norm_proj_d);
                red_param_computed = true;
            }
            a1 = std::min(a1,stepsize*delta*delta);
            aa = 2e-12*stepsize;
            stepsize = std::max(aa,a1) > min_stepsize ? std::max(aa,a1) : delta*stepsize;
        }

        if (stepsize <= min_stepsize) {
            stepsize_exit = true;
            break;
        }

        // compute a new point
        for (unsigned int i=0; i<n_ws; i++) {
            j = ind_ws[i];
            x[j] = std::max(l,std::min(v[j]+stepsize*d[j],u));
        }
        if (fast_f) {
            grad_fast();
        } else {
            grad();
        }
        f = std::inner_product(g.begin(),g.end(),x.begin(),0e0)/p_exp;
        n_f++;

    }

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Fast_atvo::grad() { // gradient of the objective function
    unsigned int k_i,k_j,len;
    double g_i,degree_i;
    const double p_exp_g = p_exp-1e0;

    g.assign(n,0e0);

    k_i = 0;
    for (unsigned int i=0; i<n-1; i++) {
        degree_i = gr->degree[i];
        while (i>gr->inst.node_index[k_i] && k_i<gr->inst.n_rows-1) {
            k_i++;
        }
        if (i != gr->inst.node_index[k_i]) {
            for (unsigned int j=i+1; j<n; j++) {
                if (x[i] != x[j]) {
                    g_i = degree_i*gr->degree[j]*sign(x[i]-x[j])*std::pow(fabs(x[i]-x[j]),p_exp_g)/gr->volume;
                    g[i] -= g_i;
                    g[j] += g_i;
                }
            }
        } else {
            k_j = 0;
            for (unsigned int j=i+1; j<n; j++) {
                if (x[i] != x[j]) {
                    len = (unsigned int) gr->inst.ind_row_col[k_i].size();
                    while (j>gr->inst.ind_row_col[k_i][k_j] && k_j<len-1) {
                        k_j++;
                    }
                    if (j != gr->inst.ind_row_col[k_i][k_j]) {
                        g_i = degree_i*gr->degree[j]*sign(x[i]-x[j])*std::pow(fabs(x[i]-x[j]),p_exp_g)/gr->volume;
                        g[i] -= g_i;
                        g[j] += g_i;
                    } else {
                        g_i = gr_inst_val_g[k_i][k_j]*sign(x[i]-x[j])*std::pow(fabs(x[i]-x[j]),p_exp_g);
                        g[i] += g_i;
                        g[j] -= g_i;
                    }
                }
            }
        }
        g[i] *= p_exp;
    }
    g[n-1] *= p_exp;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Fast_atvo::grad_fast() { // gradient of the objective function (fast computation)
    unsigned int k_i,k_j,len;
    double g_i,degree_i;
    const double p_exp_g = p_exp-1e0;

    g.assign(n,0e0);
    k_i = 0;
    for (unsigned int i=0; i<n-1; i++) {
        degree_i = gr->degree[i];
        while (i>gr->inst.node_index[k_i] && k_i<gr->inst.n_rows-1) {
            k_i++;
        }
        if (!is_in_ws[i]) {
            if (i != gr->inst.node_index[k_i]) {
                for (unsigned int j=i+1; j<n; j++) {
                    if (is_in_ws[j]) {
                        if (x[i] != x[j]) {
                            g_i = degree_i*gr->degree[j]*sign(x[i]-x[j])*std::pow(fabs(x[i]-x[j]),p_exp_g)/gr->volume;
                            g[i] -= g_i;
                            g[j] += g_i;
                        }
                        if (v[i] != v[j]) {
                            g[i] += degree_i*gr->degree[j]*sign(v[i]-v[j])*std::pow(fabs(v[i]-v[j]),p_exp_g)/gr->volume;
                        }
                    }
                }
            } else {
                k_j = 0;
                len = (unsigned int) gr->inst.ind_row_col[k_i].size();
                for (unsigned int j=i+1; j<n; j++) {
                    if (is_in_ws[j]) {
                        while (j>gr->inst.ind_row_col[k_i][k_j] && k_j<len-1) {
                            k_j++;
                        }
                        if (j != gr->inst.ind_row_col[k_i][k_j]) {
                            if (x[i] != x[j]) {
                                g_i = degree_i*gr->degree[j]*sign(x[i]-x[j])*std::pow(fabs(x[i]-x[j]),p_exp_g)/gr->volume;
                                g[i] -= g_i;
                                g[j] += g_i;
                            }
                            if (v[i] != v[j]) {
                                g[i] += degree_i*gr->degree[j]*sign(v[i]-v[j])*std::pow(fabs(v[i]-v[j]),p_exp_g)/gr->volume;
                            }
                        } else {
                            if (x[i] != x[j]) {
                                g_i = gr_inst_val_g[k_i][k_j]*sign(x[i]-x[j])*std::pow(fabs(x[i]-x[j]),p_exp_g);
                                g[i] += g_i;
                                g[j] -= g_i;
                            }
                            if (v[i] != v[j]) {
                                g[i] -= gr_inst_val_g[k_i][k_j]*sign(v[i]-v[j])*std::pow(fabs(v[i]-v[j]),p_exp_g);
                            }
                        }
                    }
                }
            }
            g[i] *= p_exp;
            g[i] += g_old[i];
        } else {
            if (i != gr->inst.node_index[k_i]) {
                for (unsigned int j=i+1; j<n; j++) {
                    if (x[i] != x[j]) {
                        g_i = degree_i*gr->degree[j]*sign(x[i]-x[j])*std::pow(fabs(x[i]-x[j]),p_exp_g)/gr->volume;
                        g[i] -= g_i;
                        g[j] += g_i;
                    }
                    if (!is_in_ws[j] && v[i]!=v[j]) {
                        g[j] -= degree_i*gr->degree[j]*sign(v[i]-v[j])*std::pow(fabs(v[i]-v[j]),p_exp_g)/gr->volume;
                    }
                }
            } else {
                k_j = 0;
                len = (unsigned int) gr->inst.ind_row_col[k_i].size();
                for (unsigned int j=i+1; j<n; j++) {
                    while (j>gr->inst.ind_row_col[k_i][k_j] && k_j<len-1) {
                        k_j++;
                    }
                    if (j != gr->inst.ind_row_col[k_i][k_j]) {
                        if (x[i] != x[j]) {
                            g_i = degree_i*gr->degree[j]*sign(x[i]-x[j])*std::pow(fabs(x[i]-x[j]),p_exp_g)/gr->volume;
                            g[i] -= g_i;
                            g[j] += g_i;
                        }
                        if (!is_in_ws[j] && v[i]!=v[j]) {
                            g[j] -= degree_i*gr->degree[j]*sign(v[i]-v[j])*std::pow(fabs(v[i]-v[j]),p_exp_g)/gr->volume;
                        }

                    } else {
                        if (x[i] != x[j]) {
                            g_i = gr_inst_val_g[k_i][k_j]*sign(x[i]-x[j])*std::pow(fabs(x[i]-x[j]),p_exp_g);
                            g[i] += g_i;
                            g[j] -= g_i;
                        }
                        if (!is_in_ws[j] && v[i]!=v[j]) {
                            g[j] += gr_inst_val_g[k_i][k_j]*sign(v[i]-v[j])*std::pow(fabs(v[i]-v[j]),p_exp_g);
                        }
                    }
                }
            }
            g[i] *= p_exp;
        }
    }
    g[n-1] *= p_exp;
    if (!is_in_ws[n-1]) {
        g[n-1] += g_old[n-1];
    }
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
inline double Fast_atvo::sign(const double& x) {
    return x > 0e0 ? 1e0 : (x < 0e0 ? -1e0 : 0e0); 
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void Fast_atvo::compute_communities(const std::vector<double>& x) {

    unsigned int h_i,h_j,ind_mod;
    double degree_cum,val_cum;
    std::vector<unsigned int> ind_sorted(n),ind_sorted_rev(n);
    std::vector<double> mod_array(n,0e0);
    
    c.assign(gr->n_original,0);
    iota(ind_sorted.begin(),ind_sorted.end(),0);
    stable_sort(ind_sorted.begin(),ind_sorted.end(),
         [&x](unsigned int i, unsigned int j) {return x[i] < x[j];});
    for (unsigned int i=0; i<n; i++) {
        ind_sorted_rev[ind_sorted[i]] = i;
    }

    for (unsigned int i=0; i<gr->inst.n_rows; i++) {
        h_i = ind_sorted_rev[gr->inst.node_index[i]];
        // check for diagonal elements in the weight matrix
        h_j = ind_sorted_rev[gr->inst.ind_row_col[i][0]];
        if (h_i < h_j) {
            mod_array[h_j] += gr->inst.val[i][0];
        } else if (h_i > h_j) {
            mod_array[h_i] += gr->inst.val[i][0];
        } else {
            mod_array[h_j] += 5e-1*gr->inst.val[i][0];
        }
        for (unsigned int j=1; j<gr->inst.ind_row_col[i].size(); j++) {
            h_j = ind_sorted_rev[gr->inst.ind_row_col[i][j]];
            if (h_i < h_j) {
                mod_array[h_j] += gr->inst.val[i][j];
            } else {
                mod_array[h_i] += gr->inst.val[i][j];
            }
        }
    }

    degree_cum = gr->degree[ind_sorted[0]];
    val_cum = mod_array[0];
    mod_array[0] = degree_cum*degree_cum/gr->volume - 2e0*val_cum;
    ind_mod = 0;
    modularity = mod_array[0];
    for (unsigned int i=1; i<n-1; i++) {
        degree_cum += gr->degree[ind_sorted[i]];
        val_cum += mod_array[i];
        mod_array[i] = degree_cum*degree_cum/gr->volume - 2e0*val_cum;
        if (mod_array[i] <= modularity) {
            ind_mod = i;
            modularity = mod_array[i];
        }
    }
    modularity *= -2e0/gr->volume;

    // compute the 0-1 vector
    for (unsigned int i=ind_mod+1; i<n; i++) {
        c[gr->ind_n_to_n_original[ind_sorted[i]]] = 1;
    }

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const std::vector<unsigned short int>& Fast_atvo::get_communities() {
    if (c.empty()) {
        c.assign(gr->n_original,0);
        compute_communities(x);
    }
    return c;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
double Fast_atvo::get_modularity() {
    if (c.empty()) {
        c.assign(gr->n_original,0);
        compute_communities(x);
    }
    return modularity;
}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
const std::vector<double>& Fast_atvo::get_x() {
    if (x_original.empty()) {
        x_original.assign(gr->n_original,0e0);
    }
    for (unsigned int i=0; i<n; i++) {
        x_original[gr->ind_n_to_n_original[i]] = x[i];
    }
    return x_original;
}
//-------------------------------------------------------------------------------------
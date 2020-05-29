// -------------------------------------------------------------------------
//
// This file is part of FAST-ATVO, which is a software for community
// detection in an undirected graph with non-negative weights.
//
// See the file 'README.txt' to know how to run the program.
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
// Last update of this file:
// May 29th, 2020
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

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <algorithm>
#include <math.h>
#include <string>
#include <string.h>
#include <time.h>
#include "graph.h"
#include "fast_atvo.h"

bool get_graph_from_file(const char*, Graph&);
bool get_x0_from_file(const char*, std::vector<double>&);
void print_usage();

int main(int argc, char *argv[]) {

    bool flag_graph,flag_x0;
    std::vector<double> x0;
    clock_t t_start,t_setup,t_solver;

    t_start = clock();

    if (argc<3 || argc>23) {
        print_usage();
        return 1;
    }

    // read graph from input file
    Graph gr;
    try {
        flag_graph = get_graph_from_file(argv[1],gr);
    }
    catch (unsigned int err) {
        std::cout << "error reading graph from '" << argv[1] << "' (line " << err << ")\n";
        return 1;
    }
    catch (char) {
        std::cout << "error: the graph volume must be positive\n";
        return 1;
    }
    if (flag_graph) {
        std::cout << "error opening '" << argv[1] << "', check if the file exists in the current directory\n";
        return 1;
    }

    // read starting point from input file
    x0.resize(gr.n_original);
    try {
        flag_x0 = get_x0_from_file(argv[2],x0);
    }
    catch (unsigned int err) {
        std::cout << "error reading starting point from '" << argv[2] << "' (line " << err << ")\n";
        return 1;
    }
    catch (char) {
        std::cout << "error: the length of the starting point must be equal to the number of nodes\n";
        return 1;
    }
    if (flag_x0) {
        std::cout << "error opening '" << argv[2] << "', check if the file exists in the current directory\n";
        return 1;
    }

    // set options (algorithm parameters + output files)
    fast_atvo_options opts;
    opts.p_exp = 14e-1;
    opts.ws_size = std::max(10,std::min(1000,int(round(3e-2*gr.n))));
    opts.out_it = 1;
    opts.lb = -1e0;
    opts.ub = 1e0;
    opts.perc_at_bounds = 1e0;
    opts.verbosity = 0;
    struct output_options {
        char *file_communities = NULL;
        char *file_modularity = NULL;
        char *file_optsol = NULL;
    } output_opts;
    if (argc > 3) {
        bool is_c_set,is_m_set,is_s_set,is_p_set,is_w_set,is_i_set,is_l_set,is_u_set,is_r_set,is_v_set;
        is_c_set = is_m_set = is_s_set = is_p_set = is_w_set = is_i_set = is_l_set = is_u_set = is_r_set = is_v_set = false;
        for (int i=3; i<argc; i=i+2) {
            if (argv[i][0]=='-' && strlen(argv[i])==2 && argc>=i+2) {
                switch (argv[i][1]) {
                case 'c':
                    output_opts.file_communities = argv[i+1];
                    is_c_set = true;
                    break;
                case 'm':
                    output_opts.file_modularity = argv[i+1];
                    is_m_set = true;
                    break;
                case 's':
                    output_opts.file_optsol = argv[i+1];
                    is_s_set = true;
                    break;
                default:
                    char *err_getnum;
                    double tmp;
                    switch (argv[i][1]) {
                    case 'p':
                        tmp = strtod(argv[i+1],&err_getnum);
                        if (*err_getnum || is_p_set) {
                            print_usage();
                            return 1;
                        }
                        if (tmp <= 1) {
                            std::cout << "error: the exponent parameter of the objective function must be greater than 1\n";
                            return 1;
                        }
                        opts.p_exp = tmp;
                        is_p_set = true;
                        break;
                    case 'w':
                        tmp = strtod(argv[i+1],&err_getnum);
                        if (*err_getnum || is_w_set) {
                            print_usage();
                            return 1;
                        }
                        if (tmp < 1) {
                            std::cout << "error: the maximum size of the working set in the optimization algorithm must be greater than or equal to 1\n";
                            return 1;
                        }
                        opts.ws_size = (unsigned int) floor(tmp);
                        is_w_set = true;
                        break;
                    case 'i':
                        tmp = strtod(argv[i+1],&err_getnum);
                        if (*err_getnum || is_i_set) {
                            print_usage();
                            return 1;
                        }
                        if (tmp < 1) {
                            std::cout << "error: the number of outer iterations for the globalization strategy must be greater than or equal to 1\n";
                            return 1;
                        }
                        opts.out_it = (unsigned int) floor(tmp);
                        is_i_set = true;
                        break;
                    case 'l':
                        tmp = strtod(argv[i+1],&err_getnum);
                        if (*err_getnum || is_l_set) {
                            print_usage();
                            return 1;
                        }
                        if (tmp >= 0) {
                            std::cout << "error: the lower bound on the variables for the optimization problem must be negative\n";
                            return 1;
                        }
                        opts.lb = tmp;
                        is_l_set = true;
                        break;
                    case 'u':
                        tmp = strtod(argv[i+1],&err_getnum);
                        if (*err_getnum || is_u_set) {
                            print_usage();
                            return 1;
                        }
                        if (tmp <= 0) {
                            std::cout << "error: the upper bound on the variables for the optimization problem must be positive\n";
                            return 1;
                        }
                        opts.ub = tmp;
                        is_u_set = true;
                        break;
                    case 'r':
                        tmp = strtod(argv[i+1],&err_getnum);
                        if (*err_getnum || is_r_set) {
                            print_usage();
                            return 1;
                        }
                        if (tmp<0e0 || tmp>1e0) {
                            std::cout << "error: the percentage of variables that will be set to the bounds in x0 must be between 0 and 1\n";
                            return 1;
                        }
                        opts.perc_at_bounds = tmp;
                        is_r_set = true;
                        break;
                    case 'v':
                        tmp = strtod(argv[i+1],&err_getnum);
                        if (*err_getnum || is_v_set) {
                            print_usage();
                            return 1;
                        }
                        if (tmp<0e0 || tmp>2e0) {
                            std::cout << "error: the verbosity level must be between 0 and 2\n";
                            return 1;
                        }
                        opts.verbosity = (unsigned short int) round(tmp);
                        is_v_set = true;
                        break;
                    default:
                        print_usage();
                        return 1;
                    }
                }
            } else {
                print_usage();
                return 1;
            }
        }
    }

    // call the solver
    t_setup = clock();
    Fast_atvo alg(&gr,x0,opts);
    alg.solve();
    t_solver = clock();

    // print final results to screen
    std::cout.precision(5);
    std::cout.setf(std::ios::fixed,std::ios::floatfield);
    std::cout << "*************************************************************"
              << "\n\nAlgorithm: FAST-ATVO"
              << "\n\ncommunity modularity = " << alg.get_modularity();
    std::cout.setf(std::ios::scientific,std::ios::floatfield);
    std::cout << "\n\nnumber of graph nodes = " << gr.n_original << " (" << gr.n_original - gr.n << " with degree zero)"
              << "\ngraph volume = " << gr.volume
              << "\nsolving time: " << (float)(t_solver-t_setup)/CLOCKS_PER_SEC << " seconds"
              << "\nset-up time: " << (float)(t_setup-t_start)/CLOCKS_PER_SEC << " seconds"
              << "\n\n*************************************************************\n\n";

    // print final results to files (if required)
    if (output_opts.file_communities != NULL) {
        const unsigned short int *c_ptr = &(alg.get_communities()[0]);
        std::cout << "writing communities in file '" << output_opts.file_communities << "'...\n";
        std::ofstream file_output;
        file_output.open(output_opts.file_communities,std::ios::trunc);
        for (unsigned int i=0; i<gr.n_original; i++) {
            file_output << *(c_ptr+i) << "\n";
        }
        file_output.close();
        std::cout << "done\n";
    }
    if (output_opts.file_modularity != NULL) {
        std::cout << "writing modularity in file '" << output_opts.file_modularity << "'...\n";
        std::ofstream file_output;
        file_output.open(output_opts.file_modularity,std::ios::trunc);
        file_output.precision(4);
        file_output.setf(std::ios::fixed,std::ios::floatfield);
        file_output << alg.get_modularity();
        file_output.close();
        std::cout << "done\n";
    }
    if (output_opts.file_optsol != NULL) {
        const double *x_ptr = &(alg.get_x()[0]);
        std::cout << "writing solution in file '" << output_opts.file_optsol << "'...\n";
        std::ofstream file_output;
        file_output.open(output_opts.file_optsol,std::ios::trunc);
        file_output.precision(6);
        file_output.setf(std::ios::scientific,std::ios::floatfield);
        for (unsigned int i=0; i<gr.n_original; i++) {
            file_output << *(x_ptr+i) << "\n";
        }
        file_output.close();
        std::cout << "done\n";
    }

    return 0;

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
bool get_graph_from_file(const char *file_name, Graph& gr) {

    //  the output is 'false' if the file is opened correctly,
    //  otherwise the output is 'true'
    //
    //  if the file is opened correctly but contains errors,
    //  - a char exception is thrown if graph volume is zero;
    //  - in case of other errors, an unsigned int exception is thrown,
    //    corresponding to the number of the line containing the error

    const char *sep_space = " ";
    const char *sep_comma = ",";
    char *token,*err_getnum;
    bool ind_i_to_read,n_rows_decreased;
    unsigned int k,q,h,b,t;
    double ind_i,ind_i_prev,ind_j,ind_j_prev,val_ij,vol,tmp_double;
    std::vector<unsigned int> idx,ind_j_vec;
    std::vector<double> weight_matrix_val;
    std::string line;
    std::ifstream data_file;

    data_file.open(file_name,std::ios_base::in);
    if (!data_file.good()) {
        data_file.close();
        return true;
    }

    // 1st loop: read file to get some information (number of nodes,
    // number of nonzero elements of the weight matrix, etc)

    // start by reading the 1st edge

    // get the first line
    std::getline(data_file,line);
    k = 1;
    if (data_file.fail() && data_file.eof()) {
        data_file.close();
        throw k;
        return false;
    }
    token = strtok(&line[0],sep_comma);
    ind_i_prev = strtod(token,&err_getnum);
    if ((*err_getnum && !isspace(*err_getnum)) || err_getnum==token || ind_i_prev<1e0 || ind_i_prev!=floor(ind_i_prev)) {
        data_file.close();
        throw k;
        return false;
    }
    token = strtok(NULL,sep_space);
    if (token == NULL) {
        data_file.close();
        throw k;
        return false;
    }
    ind_j = tmp_double = strtod(token,&err_getnum);
    if (*err_getnum || ind_j<ind_i_prev || ind_j!=floor(ind_j)) {
        data_file.close();
        throw k;
        return false;
    }
    token = strtok(NULL,sep_space);
    if (token == NULL) {
        data_file.close();
        throw k;
        return false;
    }
    val_ij = strtod(token,&err_getnum);
    if ((*err_getnum && !isspace(*err_getnum)) || err_getnum==token || val_ij<0e0) {
        data_file.close();
        throw k;
        return false;
    }

    if (val_ij > 0e0) {
        h = q = 1;
    } else {
        h = q = 0;
    }

    n_rows_decreased = false;

    token = strtok(NULL,sep_comma);

    // read the other edges

    while (1) {

        ind_i_to_read = true;
        while (token != NULL) {
            if (ind_i_to_read) {
                ind_i = strtod(token,&err_getnum);
                if ((*err_getnum && !isspace(*err_getnum)) || err_getnum==token || ind_i<ind_i_prev || ind_i != floor(ind_i)) {
                    data_file.close();
                    throw k;
                    return false;
                }
                if (ind_i > ind_i_prev) {
                    ind_i_prev = ind_i;
                    q++;
                    ind_j_prev = 0e0;
                    n_rows_decreased = false;
                } else {
                    ind_j_prev = ind_j;
                    if (n_rows_decreased) {
                        q++;
                        n_rows_decreased = false;
                    }
                }
                token = strtok(NULL,sep_space);
                ind_i_to_read = false;
            } else {
                ind_j = strtod(token,&err_getnum);
                if (*err_getnum || ind_j<ind_i || ind_j<=ind_j_prev || ind_j != floor(ind_j)) {
                    data_file.close();
                    throw k;
                    return false;
                }
                tmp_double = std::max(tmp_double,ind_j);
                token = strtok(NULL,sep_space);
                if (token == NULL) {
                    data_file.close();
                    throw k;
                    return false;
                }
                val_ij = strtod(token,&err_getnum);
                if ((*err_getnum && !isspace(*err_getnum)) || err_getnum==token || val_ij<0e0) {
                    data_file.close();
                    throw k;
                    return false;
                }
                if (val_ij > 0e0) {
                    h++;
                } else if (ind_j_prev < 1e0) {
                    q--;
                    n_rows_decreased = true;
                }
                token = strtok(NULL,sep_comma);
                ind_i_to_read = true;
            }
        }
        if (!ind_i_to_read) {
            data_file.close();
            throw k;
            return false;
        }

        //  get a new line (if any)
        if (!data_file.eof()) {
            std::getline(data_file,line);
            if (data_file.fail() && data_file.eof()) {
                k--;
            }
            token = strtok(&line[0],sep_comma);
            k++;
        } else {
            break;
        }

    }
    data_file.close();

    gr.n_original = (unsigned int) tmp_double;
    
    gr.inst.node_index.resize(q);
    gr.inst.ind_row_col.resize(q);
    gr.inst.val.resize(q);
    gr.inst.n_rows = q;

    // 2nd loop: read file again to store indices and values of the weight matrix

    // initially, 'idx' is computed such that:
    // 'idx[i]' = '2*gr.n_original' + 1 if 'i' is an isolated node;
    // 'idx[i]' = 'gr.n_original' + 1 if 'i' is a non-isolated node,
    //            but all edges of 'i' are with nodes 'j' < 'i';
    // 'idx[i]' = 'p' otherwise, where 'p' is such that there are
    //            'gr.n_original'-'p'+1 > 0 edges from 'i' and nodes 'j' >= 'i'
    b = 2*gr.n_original;
    idx.assign(gr.n_original,b+1);
    weight_matrix_val.resize(h);
    ind_j_vec.resize(h);
    data_file.open(file_name,std::ios_base::in);
    q = t = 0;
    for (unsigned int i=0; i<k; i++) {
        std::getline(data_file,line);
        token = strtok(&line[0],sep_comma);
        ind_i_to_read = true;
        while (token != NULL) {
            if (ind_i_to_read) {
                ind_i = strtod(token,&err_getnum);
                token = strtok(NULL,sep_space);
                ind_i_to_read = false;
            } else {
                ind_j = strtod(token,&err_getnum);
                token = strtok(NULL,sep_space);
                val_ij = strtod(token,&err_getnum);
                if (val_ij > 0e0) {
                    h = (unsigned int)ind_i - 1;
                    if (idx[h] > gr.n_original) {
                        if (idx[h] > b) {
                            q++;
                        }
                        idx[h] = gr.n_original + 1;
                    }
                    idx[h]--;
                    h = (unsigned int)ind_j - 1;
                    if (idx[h] > b) {
                        idx[h] = gr.n_original + 1;
                        q++;
                    }
                    weight_matrix_val[t] = val_ij;
                    ind_j_vec[t] = h;
                    t++;
                }
                token = strtok(NULL,sep_comma);
                ind_i_to_read = true;
            }
        }
    }
    gr.n = q;
    data_file.close();

    // first use 'idx' to compute 'gr.ind_n_to_n_original', 'gr.inst.node_index'
    // and to size 'gr.inst' members, then change 'idx' so that, for every
    // non-isolated node 'i' in the original graph, we have
    // 'gr.ind_n_to_n_original[idx[i]]' = 'i'
    h = gr.n_original;
    gr.ind_n_to_n_original.resize(q);
    k = q = 0;
    for (unsigned int i=0; i<h; i++) {
        if (idx[i] <= h) {
            gr.inst.node_index[k] = q;
            gr.inst.ind_row_col[k].resize(h-idx[i]+1);
            gr.inst.val[k].resize(gr.inst.ind_row_col[k].size());
            gr.ind_n_to_n_original[q] = i;
            idx[i] = q;
            q++;
            k++;
        } else if (idx[i] <= b) {
            gr.ind_n_to_n_original[q] = i;
            idx[i] = q;
            q++;
        }
    }
    
    // finally, assign values to the remaining 'gr' members
    gr.degree.assign(gr.n,0e0);
    vol = 0e0;
    q = 0;
    for (unsigned int i=0; i<gr.inst.n_rows; i++) {
        // check for a loop
        k = idx[ind_j_vec[q]];
        if (gr.inst.node_index[i] == k) {
            val_ij = weight_matrix_val[q];
            gr.inst.ind_row_col[i][0] = k;
            gr.inst.val[i][0] = val_ij;
            gr.degree[gr.inst.node_index[i]] += val_ij;
            vol += val_ij;
            t = 1;
            q++;
        } else {
            t = 0;
        }
        for (unsigned int j=t; j<gr.inst.ind_row_col[i].size(); j++) {
            val_ij = weight_matrix_val[q];
            gr.inst.ind_row_col[i][j] = idx[ind_j_vec[q]];
            gr.inst.val[i][j] = val_ij;
            gr.degree[gr.inst.node_index[i]] += val_ij;
            gr.degree[gr.inst.ind_row_col[i][j]] += val_ij;
            vol += 2e0*val_ij;
            q++;
        }
    }
    gr.volume = vol;

    if (vol <= 0) {
        throw ' ';
    }

    return false;

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
bool get_x0_from_file(const char *file_name, std::vector<double>& x0) {

    //  the output is 'false' if the file is opened correctly,
    //  otherwise the output is 'true'
    //
    //  if the file is opened correctly but contains errors:
    //  - a char exception is thrown if the length of the starting
    //    point is wrong;
    //  - in case of other errors, an unsigned int exception is thrown,
    //    corresponding to the number of the line containing the error
    
    unsigned int n = (unsigned int)x0.size();
    unsigned int k,q;
    char *token,*err_getnum;
    const char *sep_space = " ";
    std::string line;
    std::ifstream data_file;

    data_file.open(file_name,std::ios_base::in);
    if (!data_file.good()) {
        data_file.close();
        return true;
    }
    std::getline(data_file,line);
    q = 0;
    k = 1;
    if (data_file.fail() && data_file.eof()) {
        data_file.close();
        throw k;
        return false;
    }
    token = strtok(&line[0],sep_space);
    x0[q] = strtod(token,&err_getnum);
    if ((*err_getnum && !isspace(*err_getnum)) || err_getnum==token) {
        data_file.close();
        throw k;
        return false;
    }
    token = strtok(NULL," ");

    while (1) {

        while (token != NULL) {
            q++;
            if (q >= n) {
                throw ' ';
                return false;
            }
            x0[q] = strtod(token,&err_getnum);
            if ((*err_getnum && !isspace(*err_getnum)) || err_getnum==token) {
                data_file.close();
                throw k;
                return false;
            }
            token = strtok(NULL,sep_space);
        }

        if (!data_file.eof()) {
            std::getline(data_file,line);
            token = strtok(&line[0],sep_space);
            k++;
        } else {
            data_file.close();
            break;
        }

    }

    if (q < n-1) {
        throw ' ';
    }

    return false;

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
void print_usage() {
    std::cout << "\nUsage of FAST-ATVO via command prompt:\n\n"
              << "fast_atvo GraphFile X0File [options]\n\n"
              << "where '[options]' are optional input arguments that must have the\n"
              << "following form:\n\n"
              << "-c string\n"
              << "   It is the name of the file where the communities found by\n"
              << "   FAST-ATVO are printed as a 0-1 vector (if the file does not\n"
              << "   exist, then it will be created, whereas existing files with\n"
              << "   the same name will be overwritten).\n"
              << "   If not specified, by default no file is created.\n"
              << "-m string\n"
              << "   It is the name of the file where the modularity value of the\n"
              << "   communities found by FAST-ATVO is printed (if the file does not\n"
              << "   exist, then it will be created, whereas existing files with the\n"
              << "   same name will be overwritten).\n"
              << "   If not specified, by default no file is created.\n"
              << "-s string\n"
              << "   It is the name of the file where the solution found by the\n"
              << "   optimization algorithm is printed (if the file does not exist,\n"
              << "   then it will be created, whereas existing files with the same\n"
              << "   name will be overwritten).\n"
              << "   If not specified, by default no file is created.\n"
              << "-p number greater than 1\n"
              << "   It is the exponent parameter of the objective function.\n"
              << "   If not specified, by default it is equal to 1.4.\n"
              << "-w number greater than or equal to 1\n"
              << "   It is the maximum size of the working set in the optimization\n"
              << "   algorithm.\n"
              << "   If not specified, by default it is equal to\n"
              << "   max(10,min(1000,0.03*n)), where 'n' is the number of\n"
              << "   non-isolated nodes.\n"
              << "-i number greater than or equal to 1\n"
              << "   It is the number of outer iterations for the globalization\n"
              << "   strategy.\n"
              << "   If not specified, by default it is equal to 1, i.e., the\n"
              << "   globalization strategy is not activated.\n"
              << "-l number less than 0\n"
              << "   It is the lower bound on the variables for the optimization\n"
              << "   problem.\n"
              << "   If not specified, by default it is equal to -1.\n"
              << "-u number greater than 0\n"
              << "   It is the upper bound on the variables for the optimization\n"
              << "   problem.\n"
              << "   If not specified, by default it is equal to 1.\n"
              << "-r number between 0 and 1\n"
              << "   It is the percentage of negative and positive variables that\n"
              << "   will be set to the lower and upper bound in x0, respectively.\n"
              << "   If not specified, by default it is equal to 1, i.e., all\n"
              << "   non-zero variables in x0 will be set to the bounds.\n"
              << "-v number between 0 and 2\n"
              << "   It is the verbosity level, to print iteration details of the\n"
              << "   optimization algorithm.\n"
              << "   If not specified, by default it is equal to 0, i.e., there are\n"
              << "   no prints.\n\n"
              << "See the file 'README.txt' for further details.\n";
}
//-------------------------------------------------------------------------------------
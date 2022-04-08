// -------------------------------------------------------------------------
//
// This file is part of FAST-ATVO, which is a solver for community
// detection problems in undirected graphs with non-negative weights.
//
// -------------------------------------------------------------------------
//
// Reference paper:
//
// A. Cristofari, F. Rinaldi, F. Tudisco (2020). Total Variation Based
// Community Detection Using a Nonlinear Optimization Approach. SIAM Journal
// on Applied Mathematics, 80(3), 1392-1419
//
// -------------------------------------------------------------------------
//
// Authors:
// Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
// Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
// Francesco Tudisco (e-mail: francesco.tudisco@gssi.it)
//
// Last update of this file:
// April 8th, 2022
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
// Copyright 2019-2022 Andrea Cristofari, Francesco Rinaldi, Francesco
// Tudisco.
//
// -------------------------------------------------------------------------

#include "mex.h"
#include <string>
#include <algorithm>
#include <numeric>
#include <math.h>
#include "../graph.h"
#include "../fast_atvo.h"

#ifndef mxGetDoubles
#define mxGetDoubles mxGetPr
typedef double mxDouble;
#endif

// This is a MEX file for Matlab.
// See the file 'README.md' to know how to build the MEX file and run the program.

unsigned short int get_sparse_graph_matlab(Graph&, const mxArray*);
unsigned short int get_full_graph_matlab(Graph&, const mxArray*);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    unsigned short int flag_graph;
    const unsigned short int *usint_ptr;
    const double *cdbl_ptr;
    std::vector<double> x0;
    mxDouble *mdbl_ptr;

    // check the number of inupts and outputs
    if (nrhs<2) {
        mexErrMsgTxt("At least two inputs are required.");
    }
    if (nrhs>3) {
        mexErrMsgTxt("At most three inputs are required.");
    }
    if (nlhs<1) {
        mexErrMsgTxt("At least one input is required.");
    }
    if (nlhs>3) {
        mexErrMsgTxt("At most three outputs are required.");
    }
    
    if (mxIsScalar(prhs[0]) || mxGetNumberOfDimensions(prhs[0])>2 || mxGetM(prhs[0])!=mxGetN(prhs[0]) ||
        !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgTxt("The first input must be a square and upper triangular matrix with non-negative real elements.");
    }
    if (mxGetNumberOfDimensions(prhs[1])>2 || mxGetN(prhs[1])!=1 || !mxIsDouble(prhs[1]) ||
        mxIsComplex(prhs[1]) ||  mxIsSparse(prhs[1])) {
        mexErrMsgTxt("The second input must be a full real column vector with length equal to the number of nodes.");
    }
    
    // get required inputs (i.e., weight matrix and starting point)
    Graph gr;
    flag_graph = mxIsSparse(prhs[0]) ? get_sparse_graph_matlab(gr,prhs[0]) : get_full_graph_matlab(gr,prhs[0]);
    if (flag_graph > 0) {
        if (flag_graph == 1) {
            mexErrMsgTxt("The first input must be a square and upper triangular matrix with non-negative real elements.");
        } else { // i.e., flag_graph == 2
            mexErrMsgTxt("The graph volume must be positive.");
        }
    }
    if (mxGetM(prhs[1]) != gr.n_original) {
        mexErrMsgTxt("The second input must be a full real column vector of length equal to the number of nodes.");
    }
    mdbl_ptr = mxGetDoubles(prhs[1]);
    x0.resize(gr.n_original);
    for (unsigned int i=0; i<gr.n_original; i++) {
        x0[i] = mdbl_ptr[i];
    }

    // get fast-atvo options
    fast_atvo_options opts;
    opts.p_exp = 14e-1;
    opts.ws_size = std::max(10,std::min(1000,int(round(3e-2*gr.n))));
    opts.out_it = 1;
    opts.lb = -1e0;
    opts.ub = 1e0;
    opts.perc_at_bounds = 1e0;
    opts.verbosity = 0;
    if (nrhs > 2) {
        if (!mxIsStruct(prhs[2]) || mxGetNumberOfElements(prhs[2])>1) {
            mexErrMsgTxt("the third input (which is optional) must be a structure.");
        }
        for (int i=0; i<mxGetNumberOfFields(prhs[2]); i++) {
            mxArray *tmp_mxArray = mxGetFieldByNumber(prhs[2],0,i);
            const char *tmp_char = mxGetFieldNameByNumber(prhs[2],i);
            if (std::string(tmp_char).compare(std::string("p_exp")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray) || *mxGetDoubles(tmp_mxArray)<=1e0) {
                    mexErrMsgTxt("In the options, 'p_exp' must be a number greater than 1.");
                }
                opts.p_exp = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("ws_size")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray) || *mxGetDoubles(tmp_mxArray)<1e0) {
                    mexErrMsgTxt("In the options, 'ws_size' must be a number greater than or equal to 1.");
                }
                opts.ws_size = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(tmp_char).compare(std::string("out_it")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray) || *mxGetDoubles(tmp_mxArray)<1e0) {
                    mexErrMsgTxt("In the options, 'out_it' must be a number greater than or equal to 1.");
                }
                opts.out_it = (int)floor(*mxGetDoubles(tmp_mxArray));
            } else if (std::string(tmp_char).compare(std::string("lb")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray) || *mxGetDoubles(tmp_mxArray)>0e0) {
                    mexErrMsgTxt("In the options, 'lb' must be a negative number.");
                }
                opts.lb = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("ub")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray) || *mxGetDoubles(tmp_mxArray)<0e0) {
                    mexErrMsgTxt("In the options, 'ub' must be a positive number.");
                }
                opts.ub = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("perc_at_bounds")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray) ||
                    *mxGetDoubles(tmp_mxArray)<0e0 || *mxGetDoubles(tmp_mxArray)>1e0) {
                    mexErrMsgTxt("In the options, 'perc_at_bounds' must be a number between 0 and 1.");
                }
                opts.perc_at_bounds = *mxGetDoubles(tmp_mxArray);
            } else if (std::string(tmp_char).compare(std::string("verbosity")) == 0) {
                if (!mxIsScalar(tmp_mxArray) || !mxIsDouble(tmp_mxArray) || mxIsComplex(tmp_mxArray) ||
                    *mxGetDoubles(tmp_mxArray)<0e0 || *mxGetDoubles(tmp_mxArray)>2e0) {
                    mexErrMsgTxt("In the options, 'verbosity' must be a number between 0 and 2.");
                }
                opts.verbosity = (unsigned short int)round(*mxGetDoubles(tmp_mxArray));
            } else {
                mexErrMsgTxt("Not valid field name in the structure of options.");
            }            
        }
    }

    // call the solver
    Fast_atvo alg(&gr,x0,opts);
    alg.solve();

    // set outputs
    plhs[0] = mxCreateDoubleMatrix(gr.n_original,1,mxREAL);
    mdbl_ptr = mxGetDoubles(plhs[0]);
    usint_ptr = &(alg.get_communities()[0]);
    for (unsigned int i=0; i<gr.n_original; i++) {
        mdbl_ptr[i] = (double) usint_ptr[i];
    }
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleScalar(alg.get_modularity());
        if (nlhs > 2) {
            plhs[2] = mxCreateDoubleMatrix(gr.n_original,1,mxREAL);
            mdbl_ptr = mxGetDoubles(plhs[2]);
            cdbl_ptr = &(alg.get_x()[0]);
            for (unsigned int i=0; i<gr.n_original; i++) {
                mdbl_ptr[i] = cdbl_ptr[i];
            }
        }
    }
    

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
unsigned short int get_sparse_graph_matlab(Graph& gr, const mxArray *prhs) {

    // the output is:
    // 0 if no problem occurred,
    // 1 if the weight matrix is not upper triangular or has negative elements,
    // 2 if the graph volume is zero

    unsigned int k,h,q,t;
    double vol,tmp;
    double *sr;
    mwIndex *irs,*jcs;
    std::vector<unsigned int> idx,count_ind_i;
    
    h = (unsigned int) mxGetN(prhs);
    gr.n_original = h;

    sr = mxGetDoubles(prhs);
    irs = mxGetIr(prhs);
    jcs = mxGetJc(prhs);

    // initially, 'idx' is computed such that:
    // 'idx[i]' = 'gr.n_original' + 1 if 'i' is either an isolated node,
    //            or a non-isolated node linked only to nodes 'j' < 'i';
    // 'idx[i]' = 'p' otherwise, where 'p' is such that there are
    //            'gr.n_original'-'p'+1 > 0 edges from 'i' and nodes 'j' >= 'i'
    idx.assign(h,h+1);
    k = q = t = 0;
    for (unsigned int j=0; j<h; j++) {
        if (jcs[j+1]-jcs[j] > 0) {
            q++;
            for (unsigned int i=0; i<(unsigned int)(jcs[j+1]-jcs[j]); i++) {
                if (irs[t] > j) {
                    return 1;
                }
                if (idx[irs[t]] > h) {
                    k++;
                    if (jcs[irs[t]+1] == jcs[irs[t]]) {
                        q++;
                    }
                }
                idx[irs[t]]--;
                t++;
            }
        }
    }
    gr.n = q;
    gr.inst.n_rows = k;

    gr.ind_n_to_n_original.resize(q);
    gr.inst.node_index.resize(k);
    gr.inst.ind_row_col.resize(k);
    gr.inst.val.resize(k);

    // first use 'idx' to compute 'gr.ind_n_to_n_original', 'gr.inst.node_index'
    // and to size 'gr.inst' members, then change 'idx' so that, for every
    // node 'i' in the original graph linked to any node 'j' >= 'i', we have
    // 'gr.ind_n_to_n_original[gr.inst.node_index[idx[i]]]' = 'i'
    k = q = 0;
    for (unsigned int i=0; i<h; i++) {
        if (idx[i] <= h) {
            gr.inst.node_index[k] = q;
            gr.inst.ind_row_col[k].resize(h-idx[i]+1);
            gr.inst.val[k].resize(gr.inst.ind_row_col[k].size());
            gr.ind_n_to_n_original[q] = i;
            idx[i] = k;
            q++;
            k++;
        } else if (jcs[i+1]-jcs[i] > 0) {
            gr.ind_n_to_n_original[q] = i;
            q++;
        }
    }

    // finally, assign values to the remaining 'gr' members
    k = h = 0;
    count_ind_i.assign(gr.inst.n_rows,0);
    gr.degree.assign(gr.n,0e0);
    vol = 0e0;
    for (unsigned int j=0; j<gr.n_original; j++) {
        if (jcs[j+1]-jcs[j] > 0) {
            for (unsigned int i=1; i<(unsigned int)(jcs[j+1]-jcs[j]); i++) {
                if (sr[k] < 0e0) {
                    return 1;
                }
                q = idx[irs[k]];
                tmp = sr[k];
                gr.inst.ind_row_col[q][count_ind_i[q]] = h;
                gr.inst.val[q][count_ind_i[q]] = tmp;
                gr.degree[gr.inst.node_index[q]] += tmp;
                gr.degree[h] += tmp;
                vol += 2e0*tmp;
                count_ind_i[q]++;
                k++;
            }
            // check for diagonal elements in the weight matrix
            if (irs[k] == j) {
                if (sr[k] < 0e0) {
                    return 1;
                }
                q = idx[irs[k]];
                tmp = sr[k];
                gr.inst.ind_row_col[q][count_ind_i[q]] = h;
                gr.inst.val[q][count_ind_i[q]] = tmp;
                gr.degree[h] += tmp;
                vol += tmp;
            } else {
                if (irs[k]>j || sr[k]<0e0) {
                    return 1;
                }
                q = idx[irs[k]];
                tmp = sr[k];
                gr.inst.ind_row_col[q][count_ind_i[q]] = h;
                gr.inst.val[q][count_ind_i[q]] = tmp;
                gr.degree[gr.inst.node_index[q]] += tmp;
                gr.degree[h] += tmp;
                vol += 2e0*tmp;
            }
            count_ind_i[q]++;
            k++;
            h++;
        } else if (idx[j] <= gr.n_original) {
            h++;
        }
    }
    gr.volume = vol;
    
    if (vol <= 0) {
        return 2;
    }

    return 0;

}
//-------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------
unsigned short int get_full_graph_matlab(Graph& gr, const mxArray *prhs) {

    // the output is:
    // 0 if no problem occurred,
    // 1 if the weight matrix is not upper triangular or has negative elements,
    // 2 if the graph volume is zero

    unsigned int k,h,b,q;
    double vol,tmp;
    std::vector<unsigned int> idx,count_ind_i;
    mxDouble *v_ptr;

    h = (unsigned int) mxGetN(prhs);
    gr.n_original = h;

    // initially, 'idx' is computed such that:
    // 'idx[i]' = '2*gr.n_original' + 1 if 'i' is an isolated node;
    // 'idx[i]' = 'gr.n_original' + 1 if 'i' is a non-isolated node,
    //            but all edges of 'i' are with nodes 'j' < 'i';
    // 'idx[i]' = 'p' otherwise, where 'p' is such that there are
    //            'gr.n_original'-'p'+1 > 0 edges from 'i' and nodes 'j' >= 'i'
    v_ptr = mxGetDoubles(prhs);
    b = 2*h;
    idx.assign(h,b+1);
    if (*v_ptr > 0e0) {
        idx[0] = h;
        k = q = 1;
    } else if (*v_ptr < 0) {
        return 1;
    } else {
        k = q = 0;
    }
    for (unsigned int i=1; i<h; i++) {
        v_ptr++;
        if (*v_ptr != 0e0) {
            return 1;
        }
    }
    for (unsigned int j=1; j<h; j++) {
        for (unsigned int i=0; i<=j; i++) {
            v_ptr++;
            if (*v_ptr > 0e0) {
                if (idx[i] > h) {
                    if (idx[i] > b) {
                        idx[i] = h + 1;
                        q++;
                    }
                    k++;
                }
                idx[i]--;
                if (idx[j] > b) {
                    idx[j] = h + 1;
                    q++;
                }
            } else if (*v_ptr < 0) {
                return 1;
            }
        }
        for (unsigned int i=j+1; i<h; i++) {
            v_ptr++;
            if (*v_ptr != 0e0) {
                return 1;
            }
        }
    }
    gr.n = q;

    gr.inst.node_index.resize(k);
    gr.inst.ind_row_col.resize(k);
    gr.inst.val.resize(k);
    gr.inst.n_rows = k;
    gr.ind_n_to_n_original.resize(q);
    
    // first use 'idx' to compute 'gr.ind_n_to_n_original', 'gr.inst.node_index'
    // and to size 'gr.inst' members, then change 'idx' so that, for every
    // node 'i' in the original graph linked to any node 'j' >= 'i', we have
    // 'gr.ind_n_to_n_original[gr.inst.node_index[idx[i]]]' = 'i'
    k = q = 0;
    for (unsigned int i=0; i<h; i++) {
        if (idx[i] <= h) {
            gr.inst.node_index[k] = q;
            gr.inst.ind_row_col[k].resize(h-idx[i]+1);
            gr.inst.val[k].resize(gr.inst.ind_row_col[k].size());
            gr.ind_n_to_n_original[q] = i;
            idx[i] = k;
            q++;
            k++;
        } else if (idx[i] <= b) {
            gr.ind_n_to_n_original[q] = i;
            q++;
        }
    }

    // finally, assign values to the remaining 'gr' members
    h = 0;
    count_ind_i.assign(gr.inst.n_rows,0);
    v_ptr = mxGetDoubles(prhs);
    gr.degree.assign(q,0e0);
    vol = 0e0;
    if (idx[0] <= b) {
        tmp = *v_ptr;
        if (tmp > 0e0) {
            gr.inst.ind_row_col[0][0] = 0;
            gr.inst.val[0][0] = tmp;
            gr.degree[0] = tmp;
            vol = tmp;
            count_ind_i[0] = 1;
        }
        h = 1;
    }
    v_ptr += gr.n_original - 1;
    for (unsigned int j=1; j<gr.n_original; j++) {
        if (idx[j] <= b) {
            for (unsigned int i=0; i<j; i++) {
                v_ptr++;
                tmp = *v_ptr;
                if (tmp > 0e0) {
                    q = idx[i];
                    gr.inst.ind_row_col[q][count_ind_i[q]] = h;
                    gr.inst.val[q][count_ind_i[q]] = tmp;
                    gr.degree[gr.inst.node_index[q]] += tmp;
                    gr.degree[h] += tmp;
                    vol += 2e0*tmp;
                    count_ind_i[q]++;
                }
            }
            // check for diagonal elements in the weight matrix
            v_ptr++;
            tmp = *v_ptr;
            if (tmp > 0e0) {
                q = idx[j];
                gr.inst.ind_row_col[q][count_ind_i[q]] = h;
                gr.inst.val[q][count_ind_i[q]] = tmp;
                gr.degree[h] += tmp;
                vol += tmp;
                count_ind_i[q]++;
            }
            v_ptr += gr.n_original - j - 1;
            h++;
        } else {
            v_ptr += gr.n_original;
        }
    }
    gr.volume = vol;

    if (vol <= 0) {
        return 2;
    }

    return 0;

}
//-------------------------------------------------------------------------------------
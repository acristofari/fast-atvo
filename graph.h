#ifndef __GRAPH_H_INCLUDED__
#define __GRAPH_H_INCLUDED__
#include <vector>

// The weight matrix is stored as a 'weight_matrix' object named 'inst':
// - 'inst.node_index' is an array that contains all nodes 'i'
//   linked to any node 'j' >= 'i' with a positive weight
// - 'inst.ind_row_col' is an array of arrays: 'inst.ind_row_col[i]'
//   contains all nodes 'j' >= 'node_index[i]' linked to 'node_index[i]'
//   with a positive weight
// - 'inst.val' has the same structure of 'inst.ind_row_col', but
//   in each poisition ('i','j') has the weight between the corresponding
//   pair of nodes (i.e., 'inst.node_index[i]' and 'inst.ind_row_col[i][j]')
// - 'inst.n_rows' is the number of rows of the above three arrays
//
// For example, consider the following weight matrix:
//
// [0     0.9   1.5   2     0
//  0.9   0     0     0     0
//  1.5   0     0     0.8   1.1
//  2     0     0.8   0     0.3
//  0     0     1.1   0.3   0  ].
//
// Then,
//
// inst.node_index[0] = 0
// inst.node_index[1] = 2
// inst.node_index[2] = 3
//
// inst.ind_row_col[0] = [1, 2, 3]
// inst.ind_row_col[1] = [3, 4]
// inst.ind_row_col[2] = [4]
//
// inst.val[0] = [0.9, 1.5, 2]
// inst.val[1] = [0.8, 1.1]
// inst.val[2] = [0.3]
//
// inst.n_rows = 3

class Graph {
public:
    struct weight_matrix {
        std::vector<unsigned int> node_index;
        std::vector<std::vector<unsigned int>> ind_row_col;
        std::vector<std::vector<double>> val;
        unsigned int n_rows;
    } inst;
    unsigned int n,n_original; // 'n_original' is the original number of nodes,
                               // 'n' is the number of nodes with degree > 0, i.e., non-isolated nodes
                               // (nodes with degree = 0 can be ignored and a reduced graph can be considered)
    double volume;
    std::vector<unsigned int> ind_n_to_n_original; // for each node 'i' = 0,1,... in the reduced graph,
                                                   // 'ind_n_original_to_n[i]' is the index of the node
                                                   // 'i' in the original graph
    std::vector<double> degree;
};
#endif
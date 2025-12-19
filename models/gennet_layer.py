"""
PyTorch implementation of GenNet's LocallyDirected1D layer
基于稀疏连接矩阵的局部定向层，用于实现GenNet的网络结构
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
import scipy.sparse as sp
import numpy as np


class LocallyDirected1D(nn.Module):
    """
    Locally-Directed1D layer for PyTorch
    基于稀疏连接矩阵的局部定向层
    
    Args:
        mask: 稀疏连接矩阵 (scipy.sparse格式或torch.Tensor)
              shape: (n_input_nodes, n_output_nodes)
              1表示有连接，0表示无连接
        filters: 输出维度（每个输出节点的特征数）
        use_bias: 是否使用偏置
        activation: 激活函数
        kernel_initializer: 权重初始化方法
    """
    def __init__(self, mask, filters=1, use_bias=True, activation=None, 
                 kernel_initializer='xavier_uniform', dropout_rate=0.0):
        super().__init__()
        
        # 转换mask为torch格式
        if isinstance(mask, sp.spmatrix):
            # 从scipy稀疏矩阵转换为torch稀疏张量
            mask_coo = mask.tocoo()
            indices = torch.from_numpy(np.vstack([mask_coo.row, mask_coo.col])).long()
            values = torch.from_numpy(mask_coo.data).float()
            self.register_buffer('mask_sparse', torch.sparse_coo_tensor(indices, values, mask.shape))
            self.mask_dense = None
        else:
            # 如果是dense tensor
            self.register_buffer('mask_dense', mask.float())
            self.mask_sparse = None
        
        self.n_input = mask.shape[0] if isinstance(mask, sp.spmatrix) else mask.shape[0]
        self.n_output = mask.shape[1] if isinstance(mask, sp.spmatrix) else mask.shape[1]
        self.filters = filters
        self.use_bias = use_bias
        self.activation = activation
        self.dropout = nn.Dropout(dropout_rate) if dropout_rate > 0 else None
        
        # 创建权重矩阵
        # 权重形状: (n_connections, filters)
        # 其中n_connections是mask中非零元素的数量
        if isinstance(mask, sp.spmatrix):
            n_connections = mask.nnz
        else:
            n_connections = (mask > 0).sum().item()
        
        self.weight = nn.Parameter(torch.empty(n_connections, filters))
        
        if use_bias:
            self.bias = nn.Parameter(torch.zeros(self.n_output, filters))
        else:
            self.register_parameter('bias', None)
        
        # 初始化权重
        if kernel_initializer == 'xavier_uniform':
            nn.init.xavier_uniform_(self.weight)
        elif kernel_initializer == 'xavier_normal':
            nn.init.xavier_normal_(self.weight)
        elif kernel_initializer == 'kaiming_uniform':
            nn.init.kaiming_uniform_(self.weight)
        elif kernel_initializer == 'kaiming_normal':
            nn.init.kaiming_normal_(self.weight)
        else:
            nn.init.normal_(self.weight, 0, 0.01)
    
    def forward(self, x):
        """
        Forward pass
        
        Args:
            x: Input tensor, shape (batch_size, n_input, input_dim)
               或 (batch_size, n_input) 如果input_dim=1
        
        Returns:
            Output tensor, shape (batch_size, n_output, filters)
        """
        # 处理输入维度
        if x.dim() == 2:
            # (batch_size, n_input) -> (batch_size, n_input, 1)
            x = x.unsqueeze(-1)
        
        batch_size = x.size(0)
        input_dim = x.size(-1)
        
        # 获取mask
        if self.mask_sparse is not None:
            mask = self.mask_sparse.to_dense()
        else:
            mask = self.mask_dense
        
        # 更高效的实现：使用矩阵乘法
        # 方法：构建权重矩阵，然后使用矩阵乘法
        
        # 创建权重矩阵 (n_input, n_output, filters)
        weight_matrix = torch.zeros(self.n_input, self.n_output, self.filters,
                                   device=x.device, dtype=x.dtype)
        
        # 获取所有连接
        connections = torch.nonzero(mask, as_tuple=False)  # (n_connections, 2)
        
        for i, (in_idx, out_idx) in enumerate(connections):
            weight_matrix[in_idx, out_idx, :] = self.weight[i]
        
        # 如果input_dim != 1，需要先投影
        if input_dim == 1:
            x_expanded = x.squeeze(-1)  # (batch_size, n_input)
            # 使用einsum进行高效计算
            # x_expanded: (batch_size, n_input)
            # weight_matrix: (n_input, n_output, filters)
            # output: (batch_size, n_output, filters)
            output = torch.einsum('bi,iof->bof', x_expanded, weight_matrix)
        else:
            # 需要先投影input_dim到filters
            if not hasattr(self, 'input_proj'):
                # 创建投影层（如果input_dim != filters）
                self.input_proj = nn.Linear(input_dim, self.filters).to(x.device)
            
            # 投影输入: (batch_size, n_input, input_dim) -> (batch_size, n_input, filters)
            x_proj = self.input_proj(x)
            
            # 应用权重矩阵
            # x_proj: (batch_size, n_input, filters)
            # weight_matrix: (n_input, n_output, filters)
            # output: (batch_size, n_output, filters)
            output = torch.einsum('bif,iof->bof', x_proj, weight_matrix)
        
        # 添加偏置
        if self.use_bias:
            output = output + self.bias.unsqueeze(0)
        
        # 应用激活函数
        if self.activation == 'relu':
            output = F.relu(output)
        elif self.activation == 'tanh':
            output = torch.tanh(output)
        elif self.activation == 'sigmoid':
            output = torch.sigmoid(output)
        
        # 应用dropout
        if self.dropout is not None:
            output = self.dropout(output)
        
        return output


def create_topology_mask_from_file(topology_file, n_snps):
    """
    从拓扑文件创建连接矩阵
    
    Args:
        topology_file: 拓扑CSV文件路径
        n_snps: SNP数量
    
    Returns:
        scipy.sparse矩阵，shape (n_snps, n_genes)
    """
    import pandas as pd
    
    topology = pd.read_csv(topology_file)
    
    # 获取唯一的基因
    if 'layer1_node' in topology.columns:
        unique_genes = sorted(topology['layer1_node'].unique())
        n_genes = len(unique_genes)
        gene_to_idx = {gene: idx for idx, gene in enumerate(unique_genes)}
    else:
        # 单层拓扑，直接输出
        return None
    
    # 创建连接矩阵
    mask = sp.lil_matrix((n_snps, n_genes), dtype=np.float32)
    
    for _, row in topology.iterrows():
        snp_idx = int(row['layer0_node'])
        gene_idx = gene_to_idx[int(row['layer1_node'])]
        mask[snp_idx, gene_idx] = 1.0
    
    return mask.tocsr()


def create_simple_mask(n_snps, n_genes=None):
    """
    创建简单的连接矩阵（所有SNP直接连接输出）
    
    Args:
        n_snps: SNP数量
        n_genes: 基因数量（如果为None，则创建单层，直接输出）
    
    Returns:
        scipy.sparse矩阵或None
    """
    if n_genes is None:
        # 单层：所有SNP直接连接输出
        return None
    else:
        # 创建全连接矩阵（每个SNP连接到所有基因）
        mask = sp.ones((n_snps, n_genes), dtype=np.float32)
        return mask.tocsr()


"""
Multi-task Interaction Model with GenNet Integration
集成GenNet网络结构的多任务交互模型
支持两种模式：
1. 纯SNP模式：无基因拓扑，SNP直接编码
2. SNP+基因模式：使用GenNet拓扑，SNP->基因->输出
"""
import torch
import torch.nn as nn
import scipy.sparse as sp
from models.feature_selector import SparseFeatureSelector
from models.encoders_v2 import EnvironmentEncoder, MultiHeadAttention
from models.gennet_layer import LocallyDirected1D, create_topology_mask_from_file, create_simple_mask


class TaskWeightedLoss(nn.Module):
    def __init__(self, num_tasks):
        super().__init__()
        self.log_vars = nn.Parameter(torch.zeros(num_tasks))
        
    def forward(self, losses):
        weights = torch.exp(-self.log_vars)
        weighted_losses = weights * losses + self.log_vars * 0.5
        return weighted_losses.mean()


class CrossAttentionFusion(nn.Module):
    def __init__(self, hidden_dim, num_heads=8, dropout_rate=0.2):
        super().__init__()
        self.cross_attention = MultiHeadAttention(hidden_dim, num_heads)
        self.norm1 = nn.LayerNorm(hidden_dim)
        self.norm2 = nn.LayerNorm(hidden_dim)
        self.dropout = nn.Dropout(dropout_rate)
        
        self.feed_forward = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 4),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_dim * 4, hidden_dim)
        )
    
    def forward(self, geno_embed, env_embed):
        # Cross attention from genotype to environment
        geno_embed = geno_embed.unsqueeze(1)
        env_embed = env_embed.unsqueeze(1)
        
        # Genotype attends to environment
        geno_attn = self.cross_attention(geno_embed, env_embed, env_embed)
        geno_out = self.norm1(geno_embed + self.dropout(geno_attn))
        
        # Environment attends to genotype
        env_attn = self.cross_attention(env_embed, geno_embed, geno_embed)
        env_out = self.norm1(env_embed + self.dropout(env_attn))
        
        # Combine the features
        fused = geno_out + env_out
        fused = fused.squeeze(1)
        
        # Final feed-forward layer with residual connection
        ff_out = self.feed_forward(fused)
        out = self.norm2(fused + self.dropout(ff_out))
        
        return out


class TaskSpecificAttention(nn.Module):
    def __init__(self, hidden_dim, num_tasks):
        super().__init__()
        self.task_queries = nn.Parameter(torch.randn(num_tasks, hidden_dim))
        self.attention = nn.MultiheadAttention(hidden_dim, num_heads=8, batch_first=True)
        self.norm = nn.LayerNorm(hidden_dim)
        
    def forward(self, x):
        # x shape: [batch_size, hidden_dim]
        batch_size = x.size(0)
        x = x.unsqueeze(1)  # [batch_size, 1, hidden_dim]
        
        # Expand task queries for batch
        queries = self.task_queries.unsqueeze(0).expand(batch_size, -1, -1)  # [batch_size, num_tasks, hidden_dim]
        
        # Apply attention
        attn_output, _ = self.attention(queries, x, x)
        
        # Apply normalization
        attn_output = self.norm(attn_output)
        return attn_output  # [batch_size, num_tasks, hidden_dim]


class GenNetGenotypeEncoder(nn.Module):
    """
    GenNet风格的基因型编码器
    支持两种模式：
    1. 无拓扑模式：直接全连接编码
    2. 有拓扑模式：使用GenNet的LocallyDirected层
    """
    def __init__(self, geno_dim, hidden_dim, dropout_rate, 
                 topology_mask=None, use_gennet=True, num_filters=1):
        super().__init__()
        self.geno_dim = geno_dim
        self.hidden_dim = hidden_dim
        self.use_gennet = use_gennet
        self.topology_mask = topology_mask
        
        if use_gennet and topology_mask is not None:
            # 模式1: 使用GenNet拓扑（SNP -> 基因）
            self.gennet_layer = LocallyDirected1D(
                mask=topology_mask,
                filters=num_filters,
                activation='relu',
                dropout_rate=dropout_rate
            )
            
            # 获取基因数量
            if isinstance(topology_mask, sp.spmatrix):
                n_genes = topology_mask.shape[1]
            else:
                n_genes = topology_mask.shape[1]
            
            # 基因层到隐藏层的映射
            self.gene_to_hidden = nn.Sequential(
                nn.Linear(n_genes * num_filters, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout_rate)
            )
        else:
            # 模式2: 无拓扑，直接全连接编码（类似原来的GenotypeEncoder）
            self.gennet_layer = None
            self.input_layer = nn.Sequential(
                nn.Linear(geno_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout_rate)
            )
        
        # 后续的自注意力层（保持与原模型一致）
        self.self_attention = MultiHeadAttention(hidden_dim)
        self.norm1 = nn.LayerNorm(hidden_dim)
        self.norm2 = nn.LayerNorm(hidden_dim)
        
        self.feed_forward = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 4),
            nn.ReLU(),
            nn.Dropout(dropout_rate),
            nn.Linear(hidden_dim * 4, hidden_dim)
        )
    
    def forward(self, x):
        """
        Args:
            x: Input tensor, shape (batch_size, geno_dim)
        
        Returns:
            Output tensor, shape (batch_size, hidden_dim)
        """
        if self.use_gennet and self.gennet_layer is not None:
            # GenNet模式：SNP -> 基因 -> 隐藏层
            # x: (batch_size, geno_dim)
            x = x.unsqueeze(-1)  # (batch_size, geno_dim, 1) for GenNet layer
            
            # 通过GenNet层
            gene_features = self.gennet_layer(x)  # (batch_size, n_genes, filters)
            
            # 展平并映射到隐藏层
            batch_size = gene_features.size(0)
            gene_features_flat = gene_features.view(batch_size, -1)  # (batch_size, n_genes * filters)
            x = self.gene_to_hidden(gene_features_flat)  # (batch_size, hidden_dim)
        else:
            # 无拓扑模式：直接编码
            x = self.input_layer(x)  # (batch_size, hidden_dim)
        
        # 自注意力层（保持与原模型一致）
        x = x.unsqueeze(1)  # (batch_size, 1, hidden_dim)
        
        # Self attention with residual connection and layer norm
        attn_out = self.self_attention(x, x, x)
        x = self.norm1(x + attn_out)
        
        # Feed forward with residual connection and layer norm
        ff_out = self.feed_forward(x)
        x = self.norm2(x + ff_out)
        
        return x.squeeze(1)  # (batch_size, hidden_dim)


class MultiTaskGenNetModel(nn.Module):
    """
    集成GenNet的多任务交互模型
    
    Args:
        geno_dim: SNP特征维度
        env_dim: 环境特征维度
        hidden_dim: 隐藏层维度
        dropout_rate: Dropout率
        num_tasks: 任务数量
        topology_file: 拓扑文件路径（可选，如果提供则使用GenNet模式）
        use_gennet: 是否使用GenNet（如果为False，则使用简单全连接）
        num_filters: GenNet层的filter数量
    """
    def __init__(self, geno_dim, env_dim, hidden_dim, dropout_rate, num_tasks,
                 topology_file=None, use_gennet=True, num_filters=1):
        super().__init__()
        self.num_tasks = num_tasks
        self.use_gennet = use_gennet
        self.topology_file = topology_file
        
        # Feature selection
        self.feature_selector = SparseFeatureSelector(geno_dim, hidden_dim, dropout_rate)
        
        # 创建拓扑mask（如果提供拓扑文件）
        topology_mask = None
        if use_gennet and topology_file is not None:
            try:
                topology_mask = create_topology_mask_from_file(topology_file, geno_dim)
                print(f"Loaded topology from {topology_file}")
                print(f"Topology shape: {topology_mask.shape if topology_mask is not None else 'None'}")
            except Exception as e:
                print(f"Warning: Could not load topology file {topology_file}: {e}")
                print("Falling back to simple encoding mode")
                use_gennet = False
        
        # Genotype encoder (GenNet or simple)
        self.geno_encoder = GenNetGenotypeEncoder(
            geno_dim=geno_dim,
            hidden_dim=hidden_dim,
            dropout_rate=dropout_rate,
            topology_mask=topology_mask,
            use_gennet=use_gennet,
            num_filters=num_filters
        )
        
        # Environment encoder (保持不变)
        self.env_encoder = EnvironmentEncoder(env_dim, hidden_dim, dropout_rate)
        
        # Cross-attention fusion (保持不变)
        self.fusion = CrossAttentionFusion(hidden_dim, num_heads=8, dropout_rate=dropout_rate)
        
        # Task-specific attention (保持不变)
        self.task_attention = TaskSpecificAttention(hidden_dim, num_tasks)
        
        # Task-specific output layers (保持不变)
        self.task_outputs = nn.ModuleList([
            nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout_rate),
                nn.Linear(hidden_dim, 1)
            ) for _ in range(num_tasks)
        ])
        
        # Task relationships learning (保持不变)
        init_matrix = torch.eye(num_tasks) * 2.0 + torch.randn(num_tasks, num_tasks) * 0.1
        self.task_relationship = nn.Parameter(init_matrix)
        
        # Loss weighting (保持不变)
        self.task_weight_loss = TaskWeightedLoss(num_tasks)
    
    def forward(self, geno, env):
        """
        Forward pass
        
        Args:
            geno: Genotype tensor, shape (batch_size, geno_dim)
            env: Environment tensor, shape (batch_size, env_dim)
        
        Returns:
            outputs: Predictions, shape (batch_size, num_tasks)
            gates: Feature selection gates
            task_weight_loss: Task weighted loss module
        """
        # Feature selection
        geno_selected, gates = self.feature_selector(geno)
        
        # Encode genotype (GenNet or simple) and environment
        geno_embed = self.geno_encoder(geno_selected)
        env_embed = self.env_encoder(env)
        
        # Cross-attention fusion
        fused = self.fusion(geno_embed, env_embed)
        
        # Task-specific attention
        task_features = self.task_attention(fused)  # [batch_size, num_tasks, hidden_dim]
        
        # Apply task relationship learning
        task_weights = torch.softmax(self.task_relationship, dim=1)
        task_features = torch.einsum('nt,bth->bnh', task_weights, task_features)
        
        # Task-specific predictions
        outputs = []
        for i in range(self.num_tasks):
            task_output = self.task_outputs[i](task_features[:, i])
            outputs.append(task_output)
        
        outputs = torch.cat(outputs, dim=1)  # [batch_size, num_tasks]
        
        return outputs, gates, self.task_weight_loss
    
    def get_task_weights(self):
        return torch.exp(-self.task_weight_loss.log_vars)
    
    def get_task_relationships(self):
        return torch.softmax(self.task_relationship, dim=1)



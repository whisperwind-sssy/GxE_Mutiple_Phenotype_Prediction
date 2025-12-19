import torch
import torch.nn as nn
from models.feature_selector import SparseFeatureSelector
from models.encoders_v2 import GenotypeEncoder, EnvironmentEncoder, MultiHeadAttention

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

class TaskWeightedLoss(nn.Module):
    def __init__(self, num_tasks):
        super().__init__()
        self.log_vars = nn.Parameter(torch.zeros(num_tasks))
        
    def forward(self, losses):
        weights = torch.exp(-self.log_vars)
        weighted_losses = 0.5 * weights * losses + self.log_vars * 0.5
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

class MultiTaskInteractionModel(nn.Module):
    def __init__(self, geno_dim, env_dim, hidden_dim, dropout_rate, num_tasks):
        super().__init__()
        self.num_tasks = num_tasks
        
        # Feature selection and encoding
        self.feature_selector = SparseFeatureSelector(geno_dim, hidden_dim, dropout_rate)
        self.geno_encoder = GenotypeEncoder(geno_dim, hidden_dim, dropout_rate)
        self.env_encoder = EnvironmentEncoder(env_dim, hidden_dim, dropout_rate)
        
        # Cross-attention fusion
        self.fusion = CrossAttentionFusion(hidden_dim, num_heads=8, dropout_rate=dropout_rate)
        
        # Task-specific attention
        self.task_attention = TaskSpecificAttention(hidden_dim, num_tasks)
        
        # Task-specific output layers
        self.task_outputs = nn.ModuleList([
            nn.Sequential(
                nn.Linear(hidden_dim, hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout_rate),
                nn.Linear(hidden_dim, 1)
            ) for _ in range(num_tasks)
        ])
        
        # Task relationships learning
        # Initialize with small random values to encourage learning
        # This breaks symmetry and allows the model to learn task relationships
        # The initialization: identity * 2.0 gives diagonal values of 2.0
        # After softmax, this gives ~0.88 for diagonal, ~0.02 for off-diagonal
        # Adding small noise breaks symmetry and encourages learning
        init_matrix = torch.eye(num_tasks) * 2.0 + torch.randn(num_tasks, num_tasks) * 0.1
        self.task_relationship = nn.Parameter(init_matrix)
        
        # Loss weighting
        self.task_weight_loss = TaskWeightedLoss(num_tasks)

    def forward(self, geno, env):
        # Feature selection
        geno_selected, gates = self.feature_selector(geno)
        
        # Encode genotype and environment
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
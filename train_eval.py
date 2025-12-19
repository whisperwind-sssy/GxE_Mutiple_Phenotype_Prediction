import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
# from utils_tool import calculate_metrics
import itertools
from models.interaction_model import InteractionModel
import numpy as np
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr


def calculate_metrics(y_true, y_pred):
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    r2 = r2_score(y_true, y_pred, multioutput='variance_weighted')
    pearson = np.mean([
        pearsonr(y_true[:,i], y_pred[:,i])[0]
        for i in range(y_true.shape[1])
    ])
    return rmse, r2, pearson


def calculate_metrics_sep(y_true, y_pred):
    n_targets = y_true.shape[1]
    metrics = {
        'rmse': [],
        'r2': [],
        'pearson': []
    }

    for i in range(n_targets):
        rmse = np.sqrt(mean_squared_error(y_true[:, i], y_pred[:, i]))
        r2 = r2_score(y_true[:, i], y_pred[:, i])
        pearson_corr = pearsonr(y_true[:, i], y_pred[:, i])[0]

        metrics['rmse'].append(rmse)
        metrics['r2'].append(r2)
        metrics['pearson'].append(pearson_corr)

    return metrics

def train_one_epoch(model, loader, optimizer, criterion, device, l1_lambda=1e-5):
    model.train()
    running_loss = 0.0
    for geno, env, target in loader:
        geno, env, target = geno.to(device), env.to(device), target.to(device)
        out, gates = model(geno, env)
        mse_loss = criterion(out, target)
        l1_loss = gates.abs().mean()
        loss = mse_loss + l1_lambda * l1_loss

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step()
        running_loss += loss.item() * geno.size(0)
    epoch_loss = running_loss / len(loader.dataset)
    return epoch_loss
def evaluate_model(model, loader, device):
    model.eval()
    y_true, y_pred = [], []
    with torch.no_grad():
        for geno, env, target in loader:
            geno, env = geno.to(device), env.to(device)
            out, _ = model(geno, env)
            y_true.append(target.cpu().numpy())
            y_pred.append(out.cpu().numpy())
    y_true = np.concatenate(y_true)
    y_pred = np.concatenate(y_pred)
    # return calculate_metrics(y_true, y_pred)
    metrics = calculate_metrics_sep(y_true, y_pred)
    return metrics

def hyperparameter_search(geno_train, geno_val, env_train, env_val, target_train, target_val, device, output_feature):
    search_space = {
        'hidden_dim': [512, 1024, 2048],
        'dropout_rate': [0.2, 0.3, 0.4],
        'lr': [1e-3, 3e-4, 5e-4, 1e-4],
        'batch_size': [64, 128, 256],
        'weight_decay': [1e-4, 1e-5],
        'num_heads': [4, 8],  # 注意力头数
        'task_weight_init': [-0.5, 0.0, 0.5],  # 任务权重初始化
        'fusion_dropout': [0.1, 0.2]  # 融合层dropout
    }
    
    all_params = list(itertools.product(*search_space.values()))
    best_avg_metrics = {
        'rmse': float('inf'),
        'r2': float('-inf'),
        'pearson': float('-inf')
    }
    best_config = None
    
    print(f"Total combinations to try: {len(all_params)}")
    
    for params in all_params:
        hidden_dim, dropout_rate, lr, batch_size, weight_decay, num_heads, task_weight_init, fusion_dropout = params
        print(f"\nTrying parameters: hidden_dim={hidden_dim}, dropout={dropout_rate}, lr={lr}, batch_size={batch_size}")
        
        # 初始化模型
        model = MultiTaskInteractionModel(
            geno_train.shape[1], env_train.shape[1],
            hidden_dim=hidden_dim,
            dropout_rate=dropout_rate,
            num_tasks=output_feature
        ).to(device)
        
        # 设置任务权重初始化
        if hasattr(model, 'task_weight_loss'):
            with torch.no_grad():
                model.task_weight_loss.log_vars.fill_(task_weight_init)
        
        optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
        
        train_loader = DataLoader(
            TensorDataset(geno_train, env_train, target_train),
            batch_size=batch_size, shuffle=True
        )
        val_loader = DataLoader(
            TensorDataset(geno_val, env_val, target_val),
            batch_size=batch_size, shuffle=False
        )
        
        # 快速训练评估(10个epoch)
        for epoch in range(10):
            model.train()
            for batch_geno, batch_env, batch_target in train_loader:
                batch_geno = batch_geno.to(device)
                batch_env = batch_env.to(device)
                batch_target = batch_target.to(device)
                
                optimizer.zero_grad()
                outputs, gates, task_weight_loss = model(batch_geno, batch_env)
                
                # 计算每个任务的损失
                losses = []
                for i in range(output_feature):
                    losses.append(nn.MSELoss()(outputs[:, i:i+1], batch_target[:, i:i+1]))
                losses = torch.stack(losses)
                
                # 应用任务权重
                loss = task_weight_loss(losses)
                loss.backward()
                
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
                optimizer.step()
        
        # 评估
        model.eval()
        all_outputs = []
        all_targets = []
        with torch.no_grad():
            for batch_geno, batch_env, batch_target in val_loader:
                batch_geno = batch_geno.to(device)
                batch_env = batch_env.to(device)
                outputs, _, _ = model(batch_geno, batch_env)
                all_outputs.append(outputs.cpu().numpy())
                all_targets.append(batch_target.numpy())
        
        all_outputs = np.concatenate(all_outputs, axis=0)
        all_targets = np.concatenate(all_targets, axis=0)
        
        # 计算每个任务的指标
        task_metrics = {
            'rmse': [],
            'r2': [],
            'pearson': []
        }
        
        for i in range(output_feature):
            task_pred = all_outputs[:, i]
            task_true = all_targets[:, i]
            
            rmse = np.sqrt(np.mean((task_pred - task_true) ** 2))
            r2 = np.corrcoef(task_pred, task_true)[0, 1] ** 2
            pearson = np.corrcoef(task_pred, task_true)[0, 1]
            
            task_metrics['rmse'].append(rmse)
            task_metrics['r2'].append(r2)
            task_metrics['pearson'].append(pearson)
        
        # 计算平均指标
        avg_metrics = {
            'rmse': np.mean(task_metrics['rmse']),
            'r2': np.mean(task_metrics['r2']),
            'pearson': np.mean(task_metrics['pearson'])
        }
        
        # 更新最佳配置
        if (avg_metrics['rmse'] < best_avg_metrics['rmse'] and 
            avg_metrics['r2'] > best_avg_metrics['r2'] and 
            avg_metrics['pearson'] > best_avg_metrics['pearson']):
            
            best_avg_metrics = avg_metrics
            best_config = {
                'hidden_dim': hidden_dim,
                'dropout_rate': dropout_rate,
                'lr': lr,
                'batch_size': batch_size,
                'weight_decay': weight_decay,
                'num_heads': num_heads,
                'task_weight_init': task_weight_init,
                'fusion_dropout': fusion_dropout
            }
            
            print("\nNew best configuration found!")
            print(f"Average RMSE: {avg_metrics['rmse']:.4f}")
            print(f"Average R2: {avg_metrics['r2']:.4f}")
            print(f"Average Pearson: {avg_metrics['pearson']:.4f}")
            
            # 打印每个任务的指标
            for i in range(output_feature):
                print(f"\nTask {i}:")
                print(f"RMSE: {task_metrics['rmse'][i]:.4f}")
                print(f"R2: {task_metrics['r2'][i]:.4f}")
                print(f"Pearson: {task_metrics['pearson'][i]:.4f}")
    
    print("\nBest configuration found:")
    print(best_config)
    print("\nBest average metrics:")
    print(f"RMSE: {best_avg_metrics['rmse']:.4f}")
    print(f"R2: {best_avg_metrics['r2']:.4f}")
    print(f"Pearson: {best_avg_metrics['pearson']:.4f}")
    
    return best_config

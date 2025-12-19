import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import itertools
from models.interaction_model_v3 import MultiTaskInteractionModel

def hyperparameter_search_v3(geno_train, geno_val, env_train, env_val, target_train, target_val, device, num_tasks):
    """
    为多任务学习模型进行超参数搜索
    
    参数:
    - geno_train, geno_val: 训练和验证集的基因型数据
    - env_train, env_val: 训练和验证集的环境数据
    - target_train, target_val: 训练和验证集的目标变量
    - device: 计算设备(CPU/GPU)
    - num_tasks: 任务数量
    
    返回:
    - best_config: 最佳超参数配置
    """
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
            num_tasks=num_tasks
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
                for i in range(num_tasks):
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
        
        for i in range(num_tasks):
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
            for i in range(num_tasks):
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

if __name__ == "__main__":
    # 这里可以添加一个简单的测试用例
    pass 
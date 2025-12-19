"""
Main training script for GenNet-integrated Multi-task Interaction Model
集成GenNet的多任务交互模型训练脚本
"""
import pandas as pd
import torch
import torch.optim as optim
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
from models.interaction_model_gennet import MultiTaskGenNetModel
from hyperparameter_search_v3 import hyperparameter_search_v3
from train_eval import train_one_epoch, evaluate_model
from data_process.data_convert import load_or_convert_data
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import seaborn as sns
import argparse

current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

def plot_task_relationships(task_relationships, task_names, save_dir):
    os.makedirs(save_dir, exist_ok=True)
    plt.figure(figsize=(8, 6))
    sns.heatmap(task_relationships.detach().cpu().numpy(), 
                annot=True, fmt='.2f', 
                xticklabels=task_names,
                yticklabels=task_names)
    plt.title('Task Relationships')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'task_relationships_gennet_{current_time}.png'))
    plt.close()

def plot_task_weights(task_weights, task_names, save_dir):
    os.makedirs(save_dir, exist_ok=True)
    plt.figure(figsize=(10, 4))
    plt.bar(task_names, task_weights.detach().cpu().numpy())
    plt.title('Task Weights')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'task_weights_gennet_{current_time}.png'))
    plt.close()

def plot_per_trait_curves(train_losses, val_metrics_list, save_dir, trait_names):
    os.makedirs(save_dir, exist_ok=True)
    epochs = range(1, len(train_losses) + 1)

    # Plot training loss
    plt.figure(figsize=(8, 6))
    plt.plot(epochs, train_losses, label='Train Loss', color='purple')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Training Loss (GenNet Model)')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'train_loss_gennet_{current_time}.png'))
    plt.close()

    # Plot metrics for each trait
    for i, trait in enumerate(trait_names):
        plt.figure(figsize=(10, 6))
        
        rmse_curve = [epoch_vals['rmse'][i] for epoch_vals in val_metrics_list]
        r2_curve = [epoch_vals['r2'][i] for epoch_vals in val_metrics_list]
        pearson_curve = [epoch_vals['pearson'][i] for epoch_vals in val_metrics_list]

        plt.plot(epochs, rmse_curve, label='RMSE', color='r')
        plt.plot(epochs, r2_curve, label='R2', color='g')
        plt.plot(epochs, pearson_curve, label='Pearson', color='b')

        plt.xlabel('Epoch')
        plt.title(f'Metrics for {trait} (GenNet Model)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, f'trait_{trait}_metrics_gennet_{current_time}.png'))
        plt.close()

def main():
    parser = argparse.ArgumentParser(description='Train GenNet-integrated Multi-task Model')
    parser.add_argument('--topology_file', type=str, default=None,
                       help='Path to topology file (CSV format). If None, uses simple encoding mode.')
    parser.add_argument('--use_gennet', action='store_true', default=False,
                       help='Use GenNet topology (requires topology_file). If False, uses simple encoding.')
    parser.add_argument('--num_filters', type=int, default=1,
                       help='Number of filters for GenNet layer (default: 1)')
    parser.add_argument('--epochs', type=int, default=100,
                       help='Number of training epochs (default: 100)')
    parser.add_argument('--hyperparameter_search', action='store_true', default=False,
                       help='Run hyperparameter search before training')
    
    args = parser.parse_args()
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    print(f"GenNet mode: {'Enabled' if args.use_gennet else 'Disabled'}")
    if args.topology_file:
        print(f"Topology file: {args.topology_file}")
    else:
        print("Topology file: None (using simple encoding mode)")

    # Data paths
    geno_file_path_csv = "data/aligned_GENO.csv"
    geno_file_path_par = "data/aligned_GENO.parquet"
    pheno_file_path_csv = "data/South_PHENO.csv"
    pheno_file_path_par = "data/South_PHENO.parquet"
    ecov_file_path_csv = "data/aligned_ECOV.csv"
    ecov_file_path_par = "data/aligned_ECOV.parquet"

    # Load data
    env_data = load_or_convert_data(ecov_file_path_csv, ecov_file_path_par)
    print("Loaded environment data")
    snp_data = load_or_convert_data(geno_file_path_csv, geno_file_path_par)
    print("Loaded SNP data")
    y_data = load_or_convert_data(pheno_file_path_csv, pheno_file_path_par)
    print("Loaded phenotype data")

    target_columns = ['yield', 'silking', 'anthesis', 'plants_stand', 'grain_moisture', 'anthesis_GDD', 'silking_GDD']
    assert all(col in y_data.columns for col in target_columns), "Missing target columns"
    y_data = y_data[target_columns]
    num_tasks = len(target_columns)

    # 新增功能：同步删除目标变量中包含NaN的行
    print(f"Original data shape: {len(y_data)} samples")

    # 创建目标变量中的缺失值掩码 - 检查所有目标列
    missing_mask = y_data[target_columns].isna().any(axis=1)
    missing_count = missing_mask.sum()

    if missing_count > 0:
        print(f"Found {missing_count} samples with missing target values. Removing these samples...")

        # 使用相同的索引同时删除三个数据集中的对应行
        common_index = y_data.index
        y_data = y_data[~missing_mask]
        env_data = env_data.loc[common_index[~missing_mask]]
        snp_data = snp_data.loc[common_index[~missing_mask]]

        print(f"Data after cleaning: {len(y_data)} samples remain")
    else:
        print("No missing values found in target columns")

    # Data preprocessing
    env_data = env_data.drop(env_data.columns[0], axis=1)
    env_data = env_data.reset_index(drop=True)
    env_data = env_data.to_numpy()

    snp_data = snp_data.drop(snp_data.columns[0], axis=1)
    snp_data = snp_data.reset_index(drop=True)
    snp_data = snp_data.to_numpy()

    y_data = y_data.reset_index(drop=True)
    y_data = y_data.to_numpy()

    # Normalize targets
    y_mean = np.mean(y_data, axis=0)
    y_std = np.std(y_data, axis=0)
    y_data = (y_data - y_mean) / y_std

    # Convert to tensors
    geno = torch.tensor(snp_data, dtype=torch.float32)
    env = torch.tensor(env_data, dtype=torch.float32)
    target = torch.tensor(y_data, dtype=torch.float32)

    # Split data
    from sklearn.model_selection import train_test_split
    geno_trainval, geno_test, env_trainval, env_test, target_trainval, target_test = train_test_split(
        geno, env, target, test_size=0.2, random_state=42
    )
    geno_train, geno_val, env_train, env_val, target_train, target_val = train_test_split(
        geno_trainval, env_trainval, target_trainval, test_size=0.2, random_state=42
    )

    # Hyperparameter search (optional)
    if args.hyperparameter_search:
        print("Starting hyperparameter search...")
        best_config = hyperparameter_search_v3(
            geno_train, geno_val,
            env_train, env_val,
            target_train, target_val,
            device, num_tasks
        )
        print("Best Hyperparameters:", best_config)
    else:
        epochs = args.epochs

    # Model configuration
    config = {
        'hidden_dim': 2048,
        'dropout_rate': 0.2,
        'lr': 0.0003,
        'batch_size': 128,
        'weight_decay': 1e-5,
        'task_rel_lr_multiplier': 10.0  # Learning rate multiplier for task_relationship (10x default)
    }
    print("\nModel configuration:", config)
    print(f"Topology file: {args.topology_file}")
    print(f"Use GenNet: {args.use_gennet}")
    print(f"Number of filters: {args.num_filters}")

    # Initialize model with GenNet
    model = MultiTaskGenNetModel(
        geno_dim=geno.shape[1],
        env_dim=env.shape[1],
        hidden_dim=config['hidden_dim'],
        dropout_rate=config['dropout_rate'],
        num_tasks=num_tasks,
        topology_file=args.topology_file,
        use_gennet=args.use_gennet,
        num_filters=args.num_filters
    ).to(device)
    
    print(f"\nModel created successfully!")
    print(f"Total parameters: {sum(p.numel() for p in model.parameters()):,}")
    print(f"Trainable parameters: {sum(p.numel() for p in model.parameters() if p.requires_grad):,}")

    # Optimizer
    # Use higher learning rate for task_relationship to encourage learning
    # Separate parameter groups for different learning rates
    task_rel_lr_multiplier = config.get('task_rel_lr_multiplier', 10.0)  # Default 10x
    optimizer = optim.AdamW([
        {
            'params': [p for n, p in model.named_parameters() if 'task_relationship' not in n],
            'lr': config['lr'],
            'weight_decay': config['weight_decay']
        },
        {
            'params': [model.task_relationship],
            'lr': config['lr'] * task_rel_lr_multiplier,
            'weight_decay': config['weight_decay']
        }
    ])
    print(f"Task relationship learning rate: {config['lr'] * task_rel_lr_multiplier:.6f} (base lr * {task_rel_lr_multiplier})")
    
    # Learning rate scheduler
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5, patience=10
    )

    # Data loaders
    train_loader = DataLoader(
        TensorDataset(geno_trainval, env_trainval, target_trainval),
        batch_size=config['batch_size'],
        shuffle=True
    )
    test_loader = DataLoader(
        TensorDataset(geno_test, env_test, target_test),
        batch_size=config['batch_size'],
        shuffle=False
    )

    # Training setup
    train_losses = []
    val_metrics_list = []
    best_metrics = {
        'rmse': np.inf * np.ones(num_tasks),
        'r2': -np.inf * np.ones(num_tasks),
        'pearson': -np.inf * np.ones(num_tasks)
    }

    # Training loop
    for epoch in range(epochs):
        model.train()
        epoch_loss = 0
        for batch_geno, batch_env, batch_target in train_loader:
            batch_geno = batch_geno.to(device)
            batch_env = batch_env.to(device)
            batch_target = batch_target.to(device)

            optimizer.zero_grad()
            outputs, gates, task_weight_loss = model(batch_geno, batch_env)
            
            # Calculate individual losses for each task
            losses = []
            for i in range(num_tasks):
                losses.append(nn.MSELoss()(outputs[:, i:i+1], batch_target[:, i:i+1]))
            losses = torch.stack(losses)
            
            # Apply task weighting
            loss = task_weight_loss(losses)
            
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            
            epoch_loss += loss.item()

        avg_epoch_loss = epoch_loss / len(train_loader)
        train_losses.append(avg_epoch_loss)
        
        # Evaluation
        model.eval()
        all_outputs = []
        all_targets = []
        with torch.no_grad():
            for batch_geno, batch_env, batch_target in test_loader:
                batch_geno = batch_geno.to(device)
                batch_env = batch_env.to(device)
                outputs, _, _ = model(batch_geno, batch_env)
                all_outputs.append(outputs.cpu().numpy())
                all_targets.append(batch_target.numpy())

        all_outputs = np.concatenate(all_outputs, axis=0)
        all_targets = np.concatenate(all_targets, axis=0)

        # Calculate metrics for each task
        metrics = {
            'rmse': [], 'r2': [], 'pearson': []
        }
        
        for i in range(num_tasks):
            task_pred = all_outputs[:, i]
            task_true = all_targets[:, i]
            
            rmse = np.sqrt(np.mean((task_pred - task_true) ** 2))
            r2 = np.corrcoef(task_pred, task_true)[0, 1] ** 2
            pearson = np.corrcoef(task_pred, task_true)[0, 1]
            
            metrics['rmse'].append(rmse)
            metrics['r2'].append(r2)
            metrics['pearson'].append(pearson)
            
            # Update best metrics
            if rmse < best_metrics['rmse'][i]:
                best_metrics['rmse'][i] = rmse
            if r2 > best_metrics['r2'][i]:
                best_metrics['r2'][i] = r2
            if pearson > best_metrics['pearson'][i]:
                best_metrics['pearson'][i] = pearson

        val_metrics_list.append(metrics)
        
        # Print progress
        print(f"\nEpoch {epoch+1}/{epochs}")
        print(f"Training Loss: {avg_epoch_loss:.6f}")
        for i, trait in enumerate(target_columns):
            print(f"{trait}:")
            print(f"  RMSE: {metrics['rmse'][i]:.4f}")
            print(f"  R2: {metrics['r2'][i]:.4f}")
            print(f"  Pearson: {metrics['pearson'][i]:.4f}")

        # Update learning rate
        scheduler.step(avg_epoch_loss)

        # Plot task relationships and weights every 10 epochs
        if (epoch + 1) % 10 == 0:
            task_relationships = model.get_task_relationships()
            task_weights = model.get_task_weights()
            plot_task_relationships(task_relationships, target_columns, 'results')
            plot_task_weights(task_weights, target_columns, 'results')

    # Final results
    print("\n" + "="*60)
    print("Best Results (GenNet Model):")
    print("="*60)
    for i, trait in enumerate(target_columns):
        print(f"\n{trait}:")
        print(f"  Best RMSE: {best_metrics['rmse'][i]:.4f}")
        print(f"  Best R2: {best_metrics['r2'][i]:.4f}")
        print(f"  Best Pearson: {best_metrics['pearson'][i]:.4f}")

    # Plot final results
    plot_per_trait_curves(
        train_losses,
        val_metrics_list,
        save_dir='results',
        trait_names=target_columns
    )

    # Save model
    model_save_path = f"final_model_gennet_{current_time}.pth"
    torch.save({
        'model_state_dict': model.state_dict(),
        'y_mean': y_mean,
        'y_std': y_std,
        'task_relationships': model.get_task_relationships(),
        'task_weights': model.get_task_weights(),
        'config': config,
        'topology_file': args.topology_file,
        'use_gennet': args.use_gennet,
        'num_filters': args.num_filters
    }, model_save_path)
    
    print(f"\nModel saved as {model_save_path}")
    print("="*60)

if __name__ == '__main__':
    main()



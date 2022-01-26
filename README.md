# **Multi-Omics-Integration**


## **Environment**

```
<Clone Repo>
cd Multi-omics-intergration
git clone https://github.com/Jin0331/Multi-omics-intergration.git

<GPU>
conda env create --file conda_env_gpu.yaml
conda activate multiomics

<CPU>
conda env create --file conda_env_cpu.yaml
conda activate multiomics-cpu

<tensorboard>
tensorboard --logdir=tb_log --bind_all
```

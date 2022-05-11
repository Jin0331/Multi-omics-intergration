# **Multi-Omics-Integration**

## **Workflow**

![workflow](https://user-images.githubusercontent.com/42958809/167774736-f059e43b-2de6-4cae-bc5a-dc9b02e0606a.png)


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

```

```
<Subgroup Detection>
example : python src/Multi-omics-integration-subgroup.py \
    -b /home/wmbio/WORK/gitworking/Multi-omics-intergration/ \
    -c COAD \
    -e 1000

<Analysis>
example : python src/Multi-omics-integration-analysis.py \
         -b /home/wmbio/WORK/gitworking/Multi-omics-intergration/ \
         -c COAD

```

```@ wmbio.co```

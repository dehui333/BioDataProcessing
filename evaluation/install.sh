#/bin/bash 


conda install -y -c bioconda pomoxis=0.3.10
conda install -y -c conda-forge -c bioconda dgenies~=1.4.0
conda install -y -c conda-forge biopython~=1.79
conda install -y -c conda-forge selenium~=4.4.0
conda install -y -c conda-forge firefox~=105.0
conda install -y -c conda-forge geckodriver~=0.30.0

config_path=$(python -c "import dgenies;print(dgenies.__file__[:-11] +'../etc/dgenies')")
cp application.properties $config_path

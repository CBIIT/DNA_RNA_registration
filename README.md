## DNA-RNA Registration Code

This README details the installation, usage, and instructions for the DNA-RNA registration code developed at NCI/NIH. This code enables the registration of RNA and DNA datasets acquired through CellVoyager microscopes. It begins with finding the translation transform between DAPI channels and applying the same transform to the DNA channels in the DNA dataset. Raw RNA dataset plus registered DNA spot channels will be saved to the output folder. At last, the RNA dataset metadata (.mlf file) is updated to include newly added registered DNA spot channels.

### Required Packages

This code requires the following Python packages:

- `pandas`
- `scikit-image`
- `scipy`
- `matplotlib`


## 📥 Installation and Running

### Using Conda and Pip

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/CBIIT/DNA_RNA_registration.git

3. **Create a Conda Environment**:
   ```bash
   conda create --name ch_registration python=3.8
   conda activate ch_registration
   
4. **Install Packages**:
   ```bash
   pip install pandas scikit-image scipy matplotlib

5. **Run The Script:**
   ```bash
   python channel_registration.py <moving_dataset_path> <template_dataset_path> <output_folder>

Replace placeholders with actual paths:
- <moving_dataset_path> : Raw RNA dataset path.
- <template_dataset_path> : Template or RNA dataset path.
- <output_folder> : Output folder path.

### Example

```bash
python channel_registration.py V:\Users_Data\Faisal\from_adib\registered_datasets\new_test\raw\230225-EXP020623-Plate1A-2nd_20230225_175023\AssayPlate_PerkinElmer_CellCarrier-384 V:\Users_Data\Faisal\from_adib\registered_datasets\new_test\raw\230211-EXP020623-Plate1A-1st_20230211_202736\AssayPlate_PerkinElmer_CellCarrier-384 V:\Users_Data\Faisal\from_adib\registered_datasets\new_test\output



- **Note: This code assumes that DAPI images in both datasets are saved to channel 1, RNA saved as channel 4, and DNA spots as channels 2 and 3. So, please follow this routine for your imaging.**



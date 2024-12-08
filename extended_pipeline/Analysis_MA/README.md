### [README.md](file:///Users/mario/Code/CAPTURE_analysis/extended_pipeline/Analysis_MA/README.md)

This folder contains scripts for analyzing animal poses after concatenating pose estimation data from all animals and sessions. The analysis includes various conditions such as baseline, formalin pain, saline control, SNI (neuropathic pain), and sham (control).

#### Prerequisites
- MATLAB
- SPM1d toolbox: [SPM1d GitHub Repository](https://github.com/0todd0000/spm1dmatlab)

#### Scripts Overview

1. **script_01__cluster_retrieval_with_concat_predictions.m**
   - This script retrieves clusters from concatenated predictions of animal poses.
   - It processes the pose estimation data and identifies clusters based on predefined criteria.

2. **script_01__02_run_extract_cluster_vectors.m**
   - This script extracts cluster vectors from the identified clusters.
   - It prepares the data for further analysis by extracting relevant features from each cluster.

3. **script_02__02.m**
   - This script performs initial preprocessing of the pose estimation data.
   - It aligns and normalizes the data for subsequent analysis steps.

4. **script_02__03.m**
   - This script calculates various metrics from the preprocessed data.
   - It includes calculations such as velocity, acceleration, and other kinematic parameters.

5. **script_02__04.m**
   - This script performs statistical analysis on the calculated metrics.
   - It uses the SPM1d toolbox for statistical parametric mapping to identify significant differences between conditions.

6. **script_02__05_walking_analysis.m**
   - This script analyzes walking bouts based on 2D positions (x, y) of markers.
   - It calculates walking bouts, analyzes angles in an egocentric reference, and performs kinematic analysis of gait cycles.
   - It also includes plotting and visualization of the results, as well as statistical analysis using the SPM1d toolbox.

#### Analysis Workflow
1. **Data Preprocessing**
   - Align and normalize pose estimation data.
   - Extract relevant features and metrics.

2. **Cluster Analysis**
   - Retrieve clusters from concatenated predictions.
   - Extract cluster vectors for further analysis.

3. **Walking Analysis**
   - Detect walking bouts and analyze angles.
   - Perform kinematic analysis of gait cycles.
   - Plot and visualize results.

4. **Statistical Analysis**
   - Use the SPM1d toolbox for statistical parametric mapping.
   - Identify significant differences between conditions.

#### Conditions Analyzed
- Baseline
- Formalin pain
- Saline control
- SNI (neuropathic pain)
- Sham (control)

#### Notes
- Ensure the SPM1d toolbox is installed and added to the MATLAB path for statistical analysis.
- The scripts are designed to handle large datasets from multiple animals and sessions.

#### Contact
For any questions or issues, please contact the repository maintainer.

---

This README provides a comprehensive overview of the analysis scripts and their functionalities. It also highlights the prerequisites and workflow for analyzing animal poses.